# Manual de evidencia - control POSITIVO

<!-- Plantilla. Procedimiento en agent/skills/curar-control-positivo/SKILL.md. Guardar como
     agent/workspace/positivos/<triplet_id>_<drug1>_<drug2>_<evento>.md
     Notacion: {a | b} = valores permitidos; <...> = completar; (vacio) = en blanco. -->

## Identificacion

- triplet_id: T<nnn>                <!-- estable y unico; nunca reusar ids retirados -->
- drug1: <sustancia; via>           <!-- valor exacto de ref_atc -->
- drug2: <sustancia; via>           <!-- valor exacto de ref_atc -->
- event_llt: (vacio | <valor exacto de ref_llt>)
- event_pt: (vacio | <valor exacto de ref_pt>)
- event_hlt: (vacio | <valor exacto de ref_hlt>)
- event_hlgt: (vacio | <valor exacto de ref_hlgt>)   <!-- solo el nivel mas fino; el resto (vacio) -->
- fecha de curacion: <AAAA-MM-DD>
- curador: Claude (agente curar-control-positivo)
- decision: {aceptado | rechazado | pendiente}
- confidence_level: {high | moderate}

## Hipotesis

- La coadministracion pediatrica de drug1 + drug2 causa/aumenta el riesgo de <evento>
  por una IDD (no por toxicidad aislada de un farmaco): <...>

## Decision y razonamiento

- Razon causal: <...>
- Mecanismo: <...>
- Evento atribuible a la interaccion (no a un solo farmaco): <...>
- Limitaciones: <...>

## CoReporte del triplete por etapa NICHD (FAERS)

<!-- salida de agent/tools/faers_triplet_coreport.R (par + evento en el mismo caso);
     por defecto todas las etapas >= 1 -->

| nichd | triplet_coreports |
|---|---:|
| term_neonatal | <int> |
| infancy | <int> |
| toddler | <int> |
| early_childhood | <int> |
| middle_childhood | <int> |
| early_adolescence | <int> |
| late_adolescence | <int> |

## Busqueda de evidencia

<!-- cada fila con fecha de consulta y resultado textual; >= 1 fuente primaria -->

| fecha | fuente/base | query o documento | resultado |
|---|---|---|---|
| <AAAA-MM-DD> | FDA label (biomcp) | <seccion> | <...> |
| <AAAA-MM-DD> | PubMed (biomcp) | <query> | <...> |
| <AAAA-MM-DD> | DDInter (agent/tools/ddinter_lookup.R) | <par> | <...> |

## Fuentes (triple-match edad + par + evento)

<!-- primaria: evidencia edad<21 + par + evento juntos, con frase textual verbatim.
     soporte: mecanismo o contexto de clase, etiquetado como tal (no primaria). -->

| rol | tipo | cita | PMID/DOI/URL | poblacion | frase textual de respaldo |
|---|---|---|---|---|---|
| {primaria | soporte} | <...> | <...> | <...> | <...> | <...> |

## Mapeo workbook (hoja `triplets`)

- control_type: positive
- interaction_type: {pharmacokinetic | pharmacodynamic | mixed | pharmaceutical | unknown}
- mechanism: <...>
- evidence_type: <...>
- evidence_level: (vacio -> el script lo deriva del source_type | {regulatory_label | controlled_study_or_meta | observational_study | case_series | single_case_report})
- pediatric_population: <...>
- age_range: <...>
- ontogenic_modulation: {yes | no | unknown}
- higher_risk_stages: (vacio | subconjunto coma-separado de {term_neonatal, infancy, toddler, early_childhood, middle_childhood, early_adolescence, late_adolescence})
- ontogeny_evidence: (vacio | <...>)
- confidence_level: {high | moderate}
- rationale: <...>
- comments: agent/workspace/positivos/<archivo>.md

## Filas para `sources`

<!-- >= 1 fila con el mismo triplet_id; cada cita con PMID/DOI o URL estable + fecha.
     source_type es un descriptor libre; la consolidacion infiere evidence_level por
     palabra clave: "label" -> regulatory_label; "meta-analysis"/"systematic review"/
     "randomized"/"population pharmacokinetic"/"modelling" -> controlled_study_or_meta;
     "retrospective"/"cohort"/"observational"/"clinical study"/"pharmacokinetic study"
     -> observational_study; "case series" -> case_series; "case report" ->
     single_case_report; "theoretical"/"predicted" -> theoretical. -->

| triplet_id | PMID_or_DOI | URL | citation | source_type | notes |
|---|---|---|---|---|---|
| T<nnn> | (vacio \| <PMID/DOI>) | <URL estable> | <cita verbatim con frase de respaldo> | <descriptor; incluir una palabra clave de arriba para fijar evidence_level> | <...> |
