#!/usr/bin/env python3
"""Append rows to the curation workbook while preserving its dropdowns.

Why this tool exists
--------------------
The workbook `input/ddi_reference_input.xlsx` drives its cell dropdowns with
data validations whose list source lives on ANOTHER sheet (ref_atc, ref_pt,
...). Cross-sheet list validations cannot be expressed with the base OOXML
`<dataValidation>` element, so Excel stores them as an x14 extension
(`<extLst>/<x14:dataValidations>`) at the end of each worksheet.

Neither `openxlsx` (R) nor `openpyxl` (Python) round-trips that x14 extension:
loading the whole workbook and re-saving it silently drops every cross-sheet
dropdown. Any "load-all / save-all" edit is therefore lossy for this file.

This tool instead edits the workbook as a ZIP archive: it inserts the new
`<row>` elements just before `</sheetData>` in the target sheet's XML and copies
every other archive member byte-for-byte, so the x14 validations survive. The
x14 sqref ranges already span the validated rows (e.g. 2:202), so appended rows
inherit their dropdown automatically.

Usage (run from ddi_reference_set/)
-----------------------------------
    python agent/tools/append_workbook_rows.py --input rows.json [--workbook PATH] [--dry-run] [--no-backup]

`rows.json` maps each target sheet name to a list of row objects keyed by
column header. Missing columns are left blank. Example:

    {
      "triplets": [
        {"triplet_id": "T006", "control_type": "positive",
         "drug1": "methotrexate; systemic", "drug2": "trimethoprim; systemic",
         "event_pt": "Drug level increased", "interaction_type": "pharmacokinetic",
         "confidence_level": "high", "...": "..."}
      ],
      "sources": [
        {"triplet_id": "T006", "PMID_or_DOI": "PMID 2231218",
         "URL": "https://pubmed.ncbi.nlm.nih.gov/2231218/",
         "citation": "Ferrazzini G, et al. ...", "source_type": "...", "notes": "..."}
      ]
    }

Safety
------
- Writes to a temp file and verifies it BEFORE replacing the workbook:
  the x14 validation count per edited sheet must be unchanged, and the expected
  rows must be present. On any failure the workbook is left untouched.
- A timestamped backup of the workbook is made before replacement (unless
  --no-backup).
- The `triplets` sheet keeps `triplet_id` unique: appending an id that already
  exists (or repeats within the batch) is refused.

After a successful append the curated-triplet index
`agent/workspace/curated_index.tsv` is refreshed from the workbook (unless
--no-index): one row per curated triplet with its pair, event PT and primary source.
"""

import argparse
import json
import os
import re
import shutil
import sys
import zipfile
from datetime import datetime
from xml.etree import ElementTree as ET
from xml.sax.saxutils import escape

MAIN_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
REL_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"

CURATED_INDEX_PATH = "agent/workspace/curated_index.tsv"


# --- column-letter helpers -------------------------------------------------

def col_to_num(letters):
    n = 0
    for ch in letters:
        n = n * 26 + (ord(ch) - ord("A") + 1)
    return n


def num_to_col(n):
    s = ""
    while n > 0:
        n, r = divmod(n - 1, 26)
        s = chr(ord("A") + r) + s
    return s


def split_ref(ref):
    m = re.match(r"([A-Z]+)(\d+)", ref)
    return m.group(1), int(m.group(2))


# --- workbook introspection (read-only, namespace-agnostic) ----------------

def read_shared_strings(zf):
    """Return the shared-string table as a list of plain-text values."""
    if "xl/sharedStrings.xml" not in zf.namelist():
        return []
    root = ET.fromstring(zf.read("xl/sharedStrings.xml"))
    out = []
    for si in root.findall(f"{{{MAIN_NS}}}si"):
        out.append("".join(t.text or "" for t in si.iter(f"{{{MAIN_NS}}}t")))
    return out


def resolve_sheet_paths(zf):
    """Map sheet display name -> worksheet xml path via workbook.xml + rels."""
    wb = ET.fromstring(zf.read("xl/workbook.xml"))
    rels = ET.fromstring(zf.read("xl/_rels/workbook.xml.rels"))
    rid_to_target = {}
    for rel in rels:
        rid_to_target[rel.get("Id")] = rel.get("Target")
    name_to_path = {}
    for sheet in wb.iter(f"{{{MAIN_NS}}}sheet"):
        rid = sheet.get(f"{{{REL_NS}}}id")
        target = rid_to_target.get(rid, "")
        if not target.startswith("xl/"):
            target = "xl/" + target.lstrip("/")
        name_to_path[sheet.get("name")] = target
    return name_to_path


def cell_value(c, shared):
    t = c.get("t")
    if t == "s":
        v = c.find(f"{{{MAIN_NS}}}v")
        return shared[int(v.text)] if v is not None and v.text else ""
    if t == "inlineStr":
        is_el = c.find(f"{{{MAIN_NS}}}is")
        return "".join(x.text or "" for x in is_el.iter(f"{{{MAIN_NS}}}t")) if is_el is not None else ""
    v = c.find(f"{{{MAIN_NS}}}v")
    return v.text if v is not None and v.text else ""


def read_sheet_model(xml_bytes, shared):
    """Return (header_map, max_row, values_by_col) for a worksheet.

    header_map: column-name -> column-letter (from row 1)
    values_by_col: column-letter -> {row_number: value} for data rows (>1)
    """
    root = ET.fromstring(xml_bytes)
    data = root.find(f"{{{MAIN_NS}}}sheetData")
    header_map, values_by_col, max_row = {}, {}, 1
    for row in data.findall(f"{{{MAIN_NS}}}row"):
        rnum = int(row.get("r"))
        max_row = max(max_row, rnum)
        for c in row.findall(f"{{{MAIN_NS}}}c"):
            col, _ = split_ref(c.get("r"))
            val = cell_value(c, shared)
            if rnum == 1:
                if val:
                    header_map[val] = col
            else:
                values_by_col.setdefault(col, {})[rnum] = val
    return header_map, max_row, values_by_col


# --- row XML construction (inline strings; empties omitted) -----------------

def build_cell(col, rnum, value):
    return (f'<c r="{col}{rnum}" t="inlineStr"><is>'
            f'<t xml:space="preserve">{escape(str(value))}</t></is></c>')


def build_row(rnum, ordered_cells):
    cells = "".join(build_cell(col, rnum, val) for col, val in ordered_cells
                    if val is not None and str(val) != "")
    return f'<row r="{rnum}">{cells}</row>'


def x14_count(xml_text):
    return xml_text.count("x14:dataValidations")


# --- curated-triplet index -------------------------------------------------

def write_curated_index(workbook, index_path):
    """Refresh the curated-triplet index from the workbook (source of truth)."""
    with zipfile.ZipFile(workbook) as zf:
        shared = read_shared_strings(zf)
        sheet_paths = resolve_sheet_paths(zf)
        t_hdr, t_max, t_vals = read_sheet_model(zf.read(sheet_paths["triplets"]), shared)
        if "sources" in sheet_paths:
            s_hdr, s_max, s_vals = read_sheet_model(zf.read(sheet_paths["sources"]), shared)
        else:
            s_hdr, s_max, s_vals = {}, 1, {}

    def cell(hdr, vals, name, row):
        col = hdr.get(name)
        return vals.get(col, {}).get(row, "") if col else ""

    # First source (PMID/DOI, else URL) per triplet_id.
    primary = {}
    if "triplet_id" in s_hdr:
        for row in range(2, s_max + 1):
            tid = cell(s_hdr, s_vals, "triplet_id", row)
            if tid and tid not in primary:
                primary[tid] = cell(s_hdr, s_vals, "PMID_or_DOI", row) or cell(s_hdr, s_vals, "URL", row)

    cols = ["triplet_id", "control_type", "drug1", "drug2", "event_pt"]
    lines = ["\t".join(cols + ["primary_source"])]
    for row in range(2, t_max + 1):
        tid = cell(t_hdr, t_vals, "triplet_id", row)
        if not tid:
            continue
        values = [cell(t_hdr, t_vals, c, row) for c in cols]
        values.append(primary.get(tid, ""))
        lines.append("\t".join(v.replace("\t", " ").replace("\n", " ") for v in values))

    os.makedirs(os.path.dirname(index_path) or ".", exist_ok=True)
    with open(index_path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")
    return len(lines) - 1


# --- main -----------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description="Append rows to the curation workbook, preserving dropdowns.")
    ap.add_argument("--input", required=True, help="JSON: {sheet_name: [ {column: value, ...}, ... ]}")
    ap.add_argument("--workbook", default="input/ddi_reference_input.xlsx")
    ap.add_argument("--dry-run", action="store_true", help="validate and report, write nothing")
    ap.add_argument("--no-backup", action="store_true")
    ap.add_argument("--index-path", default=CURATED_INDEX_PATH,
                    help="curated-triplet index refreshed from the workbook after a successful append")
    ap.add_argument("--no-index", action="store_true", help="skip refreshing the curated-triplet index")
    args = ap.parse_args()

    with open(args.input, encoding="utf-8") as fh:
        payload = json.load(fh)

    with zipfile.ZipFile(args.workbook) as zf:
        names = zf.namelist()
        shared = read_shared_strings(zf)
        sheet_paths = resolve_sheet_paths(zf)
        members = {n: zf.read(n) for n in names}

    edited = {}   # path -> new xml text
    summary = []
    for sheet_name, rows in payload.items():
        if not rows:
            continue
        if sheet_name not in sheet_paths:
            sys.exit(f"ERROR: sheet '{sheet_name}' not found in workbook.")
        path = sheet_paths[sheet_name]
        xml_text = members[path].decode("utf-8")
        header_map, max_row, values_by_col = read_sheet_model(members[path], shared)

        # validate columns
        for r in rows:
            unknown = [k for k in r if k not in header_map]
            if unknown:
                sys.exit(f"ERROR [{sheet_name}]: unknown column(s) {unknown}. "
                         f"Valid: {sorted(header_map)}")

        # triplet_id uniqueness guard (triplets sheet only)
        if sheet_name == "triplets" and "triplet_id" in header_map:
            id_col = header_map["triplet_id"]
            existing = set(values_by_col.get(id_col, {}).values())
            batch = [r.get("triplet_id", "") for r in rows]
            for tid in batch:
                if tid in existing:
                    sys.exit(f"ERROR: triplet_id '{tid}' already present in workbook.")
            if len(set(batch)) != len(batch):
                sys.exit(f"ERROR: duplicate triplet_id within the batch: {batch}")

        # build appended rows after the current last data row
        new_xml_rows = []
        for i, r in enumerate(rows):
            rnum = max_row + 1 + i
            ordered = [(header_map[col], val) for col, val in r.items()]
            ordered.sort(key=lambda cv: col_to_num(cv[0]))
            new_xml_rows.append(build_row(rnum, ordered))
        assert xml_text.count("</sheetData>") == 1, f"{path}: unexpected sheetData markup"
        new_text = xml_text.replace("</sheetData>", "".join(new_xml_rows) + "</sheetData>")

        # per-sheet guard: x14 validations must be untouched
        if x14_count(new_text) != x14_count(xml_text):
            sys.exit(f"ERROR [{sheet_name}]: x14 validation count changed; aborting.")
        edited[path] = new_text
        summary.append(f"  {sheet_name}: +{len(rows)} row(s) (rows "
                       f"{max_row + 1}-{max_row + len(rows)})")

    if not edited:
        sys.exit("Nothing to append (empty input).")

    print("Planned append:")
    print("\n".join(summary))
    if args.dry_run:
        print("dry-run: no file written.")
        return

    # write to temp, verify, then swap in
    tmp = args.workbook + ".tmp"
    with zipfile.ZipFile(args.workbook) as zin, \
            zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as zout:
        for item in zin.infolist():
            data = zin.read(item.filename)
            if item.filename in edited:
                data = edited[item.filename].encode("utf-8")
            zout.writestr(item, data)

    # verification guard on the temp file
    try:
        with zipfile.ZipFile(tmp) as zf:
            shared2 = read_shared_strings(zf)
            for sheet_name, rows in payload.items():
                if not rows:
                    continue
                path = sheet_paths[sheet_name]
                before = members[path].decode("utf-8")
                after = zf.read(path).decode("utf-8")
                if x14_count(after) != x14_count(before):
                    raise RuntimeError(f"{sheet_name}: x14 validations lost")
                _, _, vals = read_sheet_model(zf.read(path), shared2)
                if sheet_name == "triplets" and "triplets" in payload:
                    hdr, _, _ = read_sheet_model(zf.read(path), shared2)
                    id_col = hdr["triplet_id"]
                    ids = set(vals.get(id_col, {}).values())
                    for r in rows:
                        if r["triplet_id"] not in ids:
                            raise RuntimeError(f"{r['triplet_id']} missing after write")
    except Exception as exc:
        os.remove(tmp)
        sys.exit(f"ERROR: verification failed, workbook untouched: {exc}")

    if not args.no_backup:
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        backup = f"{args.workbook}.{stamp}.bak"
        shutil.copy2(args.workbook, backup)
        print(f"backup: {backup}")
    os.replace(tmp, args.workbook)
    print(f"OK: appended and validations preserved -> {args.workbook}")

    if not args.no_index:
        n = write_curated_index(args.workbook, args.index_path)
        print(f"index: {args.index_path} ({n} triplet(s))")


if __name__ == "__main__":
    main()
