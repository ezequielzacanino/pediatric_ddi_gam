# Guia de revision humana 

Procedimiento y rubro para auditar los dossiers que el
**agente** produce al curar tripletes (`suggested/positivos/`, `suggested/negativos/`).

Evaluar:

- Criterios de aceptacion del triplete según `INCLUSION_CRITERIA.md`.

Documentar:
- Archivo de review: `reviews/agent_review.xlsx`,
  recopila la review de todos los tripletes. Lo genera (con desplegables y la
  worklist precargada) `scripts/R/00_build_review_template.R`. Tiene tres hojas:
  - `review`: una fila por triplete (veredicto + rubro de desempeno + conteos).
  - `sources`: una fila por cita (existencia/metadatos/soporte/pertinencia),
    referida por `triplet_id`. Precargada
    desde la hoja `sources` del workbook de input (`triplet_id` + `citation` +
    `pmid_or_doi`); complpetar con los desplegables de verificacion.
  - `legend`: referencia a las escalas, taxonomia de fallos y anclaje en la literatura.

## Proceso

1. Generar la planilla una sola vez (no sobrescribe una existente):
   `& 'C:\Program Files\R\R-4.4.2\bin\Rscript.exe' scripts\R\00_build_review_template.R`.
   La hoja `review` ya viene con una fila por triplete registrado en la hoja
   `triplets` del workbook de input (`triplet_id` + ruta del dossier precargados;
   la ruta queda vacia si el triplete se cargo sin su dossier, lo que marca la
   omision).
2. Para cada fila/triplete, abrir su dossier (`suggested/.../<triplet_id>_*.md`).
3. **Verificar cada fuente en la hoja `sources`**, completar
   los desplegables (`exists`/`metadata_ok`/`supports_claim`/`relevant`/
   `pediatric`).
4. Completar el **rubro de desempeño** y emitir el **veredicto** en la hoja
   `review`. El veredicto es independiente del desempeño: el curador puede aceptar un triplete con un dossier con errores (tras verificación manual de la literatura), o descartar el triplete proveniente de un dossier sin errores pero que no cumple criterios.
5. Los conteos `n_*` de la hoja `review` se pueden derivar de la hoja `sources`; usar la hoja `sources` como fuente de verdad para
   las tasas a nivel de cita.

## Justificación

La revisión combina 2 dominios:
- La evaluacion de **generaciones con
citas / atribucion** en NLP
- Los marcos de **evaluacion clinica de LLM** por
expertos. 

Las dimensiones se mapean de la siguiente manera:

- **A. Veracidad de las fuentes** -> *reference hallucination / fabricated
  citations*. Verificar que cada PMID/DOI resuelve y que los metadatos coinciden
  con la fuente real. Es el modo de fallo mas reportado en LLM aplicados a
  literatura biomedica.
- **B. Soporte de la cita (faithfulness / attribution)** -> el benchmark **ALCE**
  (Gao et al., EMNLP 2023) formaliza la *citation recall*: cada afirmacion debe
  estar **respaldada** por la fuente citada (verificar que la fuente realmente
  dice lo que el agente le atribuye).
- **C. Pertinencia / precision de las citas** -> cada fuente citada debe ser realmente relevante y debe usarse la
  **mejor evidencia disponible** (no quedarse en un case report habiendo una RS o
  cohorte). Coincide con la jerarquia `evidence_level` del proyecto.
- **D. Validez pediatrica** -> *correctness* especifica del dominio; replica el
  criterio 1 de `INCLUSION_CRITERIA.md` (evidencia < 21 años, sin extrapolacion
  adulta).
- **E. Plausibilidad mecanistica** -> *reasoning faithfulness*: el mecanismo
  declarado debe ser correcto, coherente con `interaction_type`, y a la atribucion de la **interaccion**.
- **F. Correctitud del mapeo**, **G. Completitud (completeness / omission)**,
  **H. Calibracion** -> ejes del marco de **evaluacion humana de LLM en salud**
  (Tam et al., *npj Digital Medicine* 2024; revision en *BMC Med Inform Decis
  Mak* 2025): *correctness*, *completeness*, y la deteccion de *over-/under-
  calling* (sobre- o sub-llamado), que son las dos caras del sesgo de
  confirmacion. La calibracion verifica que `evidence_level`/`confidence_level`
  no sobre-estimen la evidencia.

Referencias:
- Gao T, Yen H, Yu J, Chen D. Enabling Large Language Models to Generate Text
  with Citations (ALCE). EMNLP 2023. arXiv:2305.14627.
- Tam TYC, et al. A framework for human evaluation of large language models in
  healthcare derived from literature review. npj Digit Med. 2024;7:258.
  doi:10.1038/s41746-024-01258-7.
- Sina Shool, et al. A systematic review of large language model (LLM) evaluations in clinical medicine. BMC Med Inform
  Decis Mak. 2025;25:117. doi:10.1186/s12911-025-02954-4.

## Escala

Cada item del rubro: **si / parcial / no** (o `na` si no aplica; desplegables en
la hoja `review`). "parcial" se reserva para cumplimiento incompleto verificable
(p.ej. 2 de 3 claims soportados). Justificar en `comments` cuando el valor no es
"si". Los conteos de fuentes (`n_fabricated`, `n_unsupported`, `n_irrelevant`...)
son enteros y permiten calcular tasas a nivel de cita, no solo de triplete.

## Negativos

Para un control negativo (`suggested/negativos/`), el objeto de la evidencia es la
**ausencia documentada de interaccion**, no su presencia. Adaptaciones (misma
hoja `review`, reinterpretando dos columnas):

- Hoja `sources`: verificar que las fuentes de **ausencia** existen y que la
  consulta citada es real y pertinente al
  par+evento.
- `pediatric_valid`: pasa a significar **plausibilidad de coadministracion**
  pediatrica del par.
- `mechanism_sound` + `interaction_attributed`: se reinterpretan como **riesgo de
  misclasificacion** (*harm/safety*): que el evento no sea ADR atribuible a un
  solo farmaco y que no sea "ausencia de evidencia por par poco estudiado"
  (Hauben 2016). Usar el tag `misclassified_negative` (o `single_drug_attribution`).

## Agregacion

Utilizar `reviews/agent_review.xlsx`; leer las hojas `review` y `sources` para el desempeño del agente sobre el set:

- **Tasas por triplete** (hoja `review`): % `verdict == accepted_as_is`,
  % `sources_all_real == no`, % `all_claims_supported == yes`,
  % `pediatric_valid == yes`, % `mechanism_sound == yes`, etc.
- **Tasas por cita** (hoja `sources`): `mean(exists == "no")` (tasa de
  fabricacion), `mean(supports_claim != "yes")` (no-soporte), `mean(relevant ==
  "no")` (no-pertinencia). Son el analogo de la *citation precision/recall* de
  ALCE para todo el set.
- **Perfil de fallos**: frecuencia de cada token de `failure_tags`, para ver
  donde falla mas el agente (p.ej. fabricacion vs sub-uso de evidencia fuerte).

Mantener `triplet_id` estable (igual que en el workbook de input) para poder
cruzar el desempeño con `evidence_level`, par o evento.
