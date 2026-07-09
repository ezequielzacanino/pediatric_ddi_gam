# Guía de curación para el humano (validación dual)

Procedimiento paso a paso para el **curador humano** del ejercicio de validación
dual. Es la **contraparte humana del skill `curar-control-positivo`**
(`../suggested/positivos/SKILL.md`): ambos curadores —humano y agente— deben correr
**el mismo SOP** para que el acuerdo medido refleje diferencias de juicio y no de
procedimiento. Los criterios de aceptación (el *qué*) viven en
[`../INCLUSION_CRITERIA.md`](../INCLUSION_CRITERIA.md) y no se repiten acá: esta guía
solo describe el *cómo* para este ejercicio.

## Reglas del ejercicio

- **Solo controles positivos.** Acá únicamente se buscan tripletes donde el evento
  es atribuible a la **interacción** del par (no a un solo fármaco). No hay
  `control_type` ni negativos.
- **Ciego.** No mires el workbook del otro curador (`dual_curation_agent.xlsx`) ni
  el set de referencia curado del proyecto padre mientras trabajás.
- **Pares asignados.** Tenés 50 pares droga-droga fijos en `dual_curation_human.xlsx`.
  La hoja `pairs` los lista; además hay **una hoja por par** (`pair_01` … `pair_50`).
- **Un par puede dar cero, uno o varios tripletes.** No asumas de antemano cuántos
  pares interactúan ni en qué proporción: trabajá **cada par por su propia
  evidencia**. Tanto un `no_triplet_found` honesto como uno o más tripletes bien
  documentados son resultados válidos; lo que se mide es tu juicio caso por caso.

## Dónde buscar evidencia

Las mismas fuentes que usa el agente, en este orden de prioridad:

1. **Etiquetas regulatorias:** FDA (DailyMed / label) y EMA (SmPC).
2. **Literatura:** PubMed / Europe PMC.
3. **Compendio / base de IDD** (si tenés acceso): BNF, Micromedex, DDInter,
   Drugs.com.

Priorizá **evidencia pediátrica directa del par** (< 21 años). Registrá fecha,
fuente, query y resultado de cada búsqueda (en tu propia bitácora; en el workbook
no hace falta hoja de fuentes en este ejercicio).

## Procedimiento por par (ejecutar en orden)

Para **cada uno de los 20 pares**:

### Paso 1 — Encuadre
Abrí la hoja del par (`pair_NN`). En las filas 1-2 ves `drug1` y `drug2` (fijos: no
los ingresás). Formulá: "¿la coadministración pediátrica de `drug1` + `drug2` causa
o aumenta el riesgo de algún evento por una IDD?". Anotá los eventos candidatos
plausibles y su mecanismo posible (`pharmacokinetic`, `pharmacodynamic`, `mixed`,
`pharmaceutical`, `unknown`).

### Paso 2 — Búsqueda de evidencia (irreductible)
Buscá en FDA/EMA + PubMed/Europe PMC (y un compendio de IDD si tenés). Recogé:
- ≥ 1 fuente que documente la **coadministración del par** (no cada fármaco por
  separado) en < 21 años.
- El mecanismo y la plausibilidad temporal/mecanística.
- Estudios grandes (RS, meta, ECA, cohorte, poblacional/PK) si existen; los
  reportes de caso valen como soporte salvo que sean la única evidencia pediátrica
  con dechallenge/rechallenge y mecanismo claro.

### Paso 3 — Filtro de criterios (todos obligatorios)
Rechazá el triplete si falla **cualquiera** (ver `../INCLUSION_CRITERIA.md`):
1. Población < 21 años con confirmación directa del par (no extrapolación adulta).
2. Coadministración de **ambos** fármacos documentada.
3. Evento atribuible a la **interacción** (no a un solo fármaco).
4. Mecanismo declarado (uno de los 5 tipos; nunca `none` en un positivo).
5. Evento mapeable a MedDRA (existe en los desplegables `event_pt`/`event_llt`).
6. ≥ 1 fuente trazable (PMID/DOI o URL estable).

### Paso 4 — Verificación FAERS por etapa (si está disponible)
Si tenés el coReporte por etapa NICHD del par, regla por defecto: aceptar solo si
todas las etapas tienen `coadmin_reports ≥ 1`, salvo decisión explícita. Documentá.

### Paso 5 — Mapeo del evento
Cargá el evento en el **nivel MedDRA más fino disponible** eligiéndolo del
desplegable de la columna correspondiente (`event_llt` / `event_pt` / `event_hlt` /
`event_hlgt`). Completá **solo** ese nivel; el script 02 lo resuelve a PT para
comparar. Los fármacos no se ingresan (son el contexto del par).

### Paso 6 — Nivel de evidencia y confianza
- `evidence_level` (controlado, de mayor a menor): `regulatory_label` >
  `controlled_study_or_meta` > `observational_study` > `case_series` >
  `single_case_report` > `theoretical`.
- `confidence_level`: `high` (etiqueta/ECA/meta/RS, o caso/serie pediátrica con
  relación temporal+mecanística clara e idealmente confirmación independiente) o
  `moderate` (caso único, mecanismo incierto, atribución menos directa). No aceptes
  positivos que solo alcanzarían nivel "bajo" (teórico puro).

### Paso 7 — Modulación ontogénica (opcional)
Si la evidencia indica que el riesgo **varía con la edad**: `ontogenic_modulation =
yes` y completá `higher_risk_stages` (subconjunto de las 7 etapas NICHD) y
`ontogeny_evidence`. Si no lo evaluaste, dejá `unknown`.

### Paso 8 — Registrar en la hoja del par
- **Encontraste uno o más tripletes:** una **fila por triplete** debajo del
  encabezado. Completá el evento (Paso 5), `interaction_type`, `mechanism`,
  `evidence_level`, `confidence_level`, las columnas de ontogenia (Paso 7) y, si
  querés dejar la cita a mano, `source_title`/`source_year`/`source_type` y
  `rationale`/`comments`. Asignale un `triplet_id` propio (p.ej. `H001`, `H002`…)
  como referencia tuya; no es necesario que coincida con el del agente (el matching
  se hace por par + evento, no por id).
- **No encontraste ninguno:** marcá `no_triplet_found = yes` en la primera fila de
  la hoja y dejá el evento en blanco.

## Exclusiones (rechazo automático)

- Evidencia solo adulta sin confirmación pediátrica del par.
- Interacción solo teórica o predicha por una base de datos sin caso real.
- Evento atribuible a un único fármaco.
- Evento no mapeable a ningún nivel MedDRA del desplegable.

## Cuando termines

Guardá `dual_curation_human.xlsx`. La comparación con el agente la corre el
investigador con `scripts/R/02_compare_dual_curation.R`; vos solo entregás el
workbook completo.

## Checklist por par

Criterios 1-6 OK → mecanismo declarado → evento en el nivel MedDRA más fino →
FAERS por etapa (si hay) → fila(s) en la hoja del par, **o** `no_triplet_found =
yes` si el par no dio ninguno.
