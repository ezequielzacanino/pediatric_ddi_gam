# Agente de curación — set de referencia pediátrico de IDD

Carpeta autocontenida para curar controles (positivos y negativos) del set de
referencia de interacciones farmaco-farmaco (IDD) en pediatría. Todo el contexto
que el agente necesita vive acá; lo que produce se guarda acá.

## Contrato de trabajo

- **Working directory:** ejecutar siempre desde la raíz del proyecto
  (`ddi_reference_set/`). Las tools resuelven sus rutas contra ese directorio,
  no contra su propia ubicación.
- **Precondición:** el workbook `input/ddi_reference_input.xlsx` ya existe (lo crea
  el humano con `scripts/R/00_build_input_template.R`). El agente lo lee y le
  agrega filas; no lo genera.

## Qué leer (contexto)

- `agent/skills/curar-control-positivo/SKILL.md` — procedimiento paso a paso (positivos).
- `agent/skills/curar-control-negativo/SKILL.md` — procedimiento paso a paso (negativos).
- `agent/skills/*/TEMPLATE.md` — plantilla del dossier de evidencia.
- `../INCLUSION_CRITERIA.md` (raíz) — criterios de aceptación (el *qué*).
- `../CURATION_WORKFLOW.md` (raíz) — flujo general (el *cómo*).

## Herramientas (correr desde `ddi_reference_set/`)

- `agent/tools/vocab_lookup.R` — resolver strings exactos de los desplegables
  (`ref_atc`, `ref_pt/llt/hlt/hlgt`) sin cargar las hojas enteras.
- `agent/tools/ddinter_lookup.R` — soporte mecanístico DDInter (CSVs en `input/ddinter/`).
- `agent/tools/faers_triplet_coreport.R` — detectabilidad del triplete (coReporte del triplete
  por etapa NICHD: par + evento en el mismo caso). Con `--mode rank` (solo los dos farmacos):
  tabla rankeada de los eventos PT coReportados con el par, con `triplet_coreports` y
  `single_drug_event_max` por evento.
- `agent/tools/append_workbook_rows.py` — agregar filas al workbook (`triplets`/`sources`)
  a partir de un JSON, preservando los desplegables. Edita el XML del `.xlsx` como ZIP:
  openxlsx (R) y openpyxl (Python) **descartan** las validaciones de lista cross-hoja
  (extensión x14) al reguardar el archivo entero, así que el volcado normal borra los
  dropdowns. El tool inserta las filas sin tocar el resto, hace backup y verifica que las
  validaciones sigan presentes (si no, no escribe). Uso:
  `python agent/tools/append_workbook_rows.py --input filas.json [--dry-run]`.
  Tras un append exitoso refresca `agent/workspace/curated_index.tsv` desde el workbook.

Evidencia primaria vía MCP (`biomcp`: PubMed, etiquetas/FAERS openFDA, ensayos),
según el gate de citas del SKILL.

## Dónde escribe (salida)

- Dossiers de evidencia: `agent/workspace/positivos/<triplet_id>_*.md` y
  `agent/workspace/negativos/<triplet_id>_*.md`.
- Filas en el workbook `input/ddi_reference_input.xlsx` (hojas `triplets` y `sources`),
  agregadas con `agent/tools/append_workbook_rows.py` para no perder los desplegables.
- Índice de tripletes curados: `agent/workspace/curated_index.tsv`, derivado del workbook
  por `agent/tools/append_workbook_rows.py` (par + evento + fuente ya presentes).

## Qué NO tocar

- `scripts/R/` — pipeline del humano (build de templates, generación de candidatos,
  consolidación del set final, review). El agente no lo ejecuta ni lo modifica.
- `results/`, `reviews/`, y el resto de subproyectos del workspace.
