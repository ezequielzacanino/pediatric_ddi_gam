# Manual de evidencia - control NEGATIVO

<!-- Plantilla. Procedimiento en agent/skills/curar-control-negativo/SKILL.md. Guardar como
     agent/workspace/negativos/<triplet_id>_<drug1>_<drug2>_<evento>.md
     Notacion: {a | b} = valores permitidos; <...> = completar; (vacio) = en blanco. -->

## Identificacion

- triplet_id: N<nnn>                <!-- estable y propio; NO el posicional del CSV -->
- drug1: <sustancia; via>           <!-- valor exacto de ref_atc -->
- drug2: <sustancia; via>           <!-- valor exacto de ref_atc -->
- event_llt: (vacio | <valor exacto de ref_llt>)
- event_pt: (vacio | <valor exacto de ref_pt>)
- event_hlt: (vacio | <valor exacto de ref_hlt>)
- event_hlgt: (vacio | <valor exacto de ref_hlgt>)   <!-- solo el nivel mas fino; el resto (vacio) -->
- fecha de curacion: <AAAA-MM-DD>
- curador: Claude (agente curar-control-negativo)
- decision: {aceptado | rechazado | pendiente}
- confidence_level: {high | moderate}

## Candidato de origen

- csv_row (cursor del harness): <int>
- candidato_posicional_CSV: N<nnn>          <!-- solo referencia; se renumera por corrida -->
- match_strategy: {event_swap | drug_swap}
- matched_positive: T<nnn>
- known_interacting_pair: {TRUE | FALSE}
- pair_coreport: <int>
- triplet_coreport: <int>

## Decision y razonamiento

- Par plausible en pediatria: <...>
- Sin interaccion documentada para este evento: <...>
- Evento no explicado por un solo farmaco (por etiqueta): <...>
- Riesgo de misclasificacion: <...>
- Limitaciones: <...>

## Busqueda de ausencia de interaccion

<!-- >= 2 fuentes independientes, cada una con fecha de consulta y resultado textual -->

| fecha | fuente/base | query o documento | resultado |
|---|---|---|---|
| <AAAA-MM-DD> | DDInter (agent/tools/ddinter_lookup.R) | <par> | <...> |
| <AAAA-MM-DD> | FDA label drug1 (biomcp) | <seccion> | <...> |
| <AAAA-MM-DD> | FDA label drug2 (biomcp) | <seccion> | <...> |
| <AAAA-MM-DD> | PubMed (biomcp) | <query> | <...> |

## Revision de ADR mono-farmaco (por etiqueta)

| farmaco | fuente | link | hallazgo sobre el evento |
|---|---|---|---|
| drug1 | <fuente> | <URL estable> | <...> |
| drug2 | <fuente> | <URL estable> | <...> |

## Mapeo workbook (hoja `triplets`)

- control_type: negative
- interaction_type: none
- mechanism: (vacio)
- ontogenic_modulation: {no | unknown}
- higher_risk_stages: (vacio)
- ontogeny_evidence: (vacio)
- evidence_type: (vacio | <valor>)
- evidence_level: (vacio)           <!-- lo deriva el script de consolidacion -->
- pediatric_population: <...>
- age_range: <...>
- source_title: <texto; menciona ausencia + fecha de consulta>
- source_year: <AAAA>
- source_type: compendia_and_labels
- rationale: <...>
- comments: agent/workspace/negativos/<archivo>.md

## Filas para `sources`

<!-- >= 1 fila; la cita menciona la ausencia y la fecha de consulta -->

| triplet_id | PMID_or_DOI | URL | citation | source_type | notes |
|---|---|---|---|---|---|
| N<nnn> | (vacio \| <PMID/DOI>) | <URL estable> | <cita: ausencia + fecha> | {interaction_compendium \| product_label \| literature \| reference} | {absence documented for this event \| single-drug attribution excluded \| no literature for pair-event association} |
