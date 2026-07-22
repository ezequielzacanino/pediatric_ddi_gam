---
name: curar-control-positivo
description: >
  Curar un CONTROL POSITIVO del set de referencia pediatrico de IDD: un triplete
  farmaco1 x farmaco2 x evento MedDRA en el que el evento es atribuible a la
  INTERACCION del par (no a un solo farmaco) en poblacion < 21 anos. Usar cuando
  se quiera proponer, evaluar o documentar un nuevo positivo, o verificar que un
  candidato cumple los criterios de inclusion antes de cargarlo al workbook.
---

# Skill: curar un control POSITIVO

Procedimiento autonomo y reproducible para convertir una hipotesis de interaccion
en un triplete positivo trazable, listo para cargar en
`input/ddi_reference_input.xlsx`. El *que* (criterios de aceptacion) vive en
`INCLUSION_CRITERIA.md`; el *como* general en `CURATION_WORKFLOW.md`. Este skill
es la version operativa paso a paso, pensada para ejecutarse fila por fila.

## Cuando aplica

Un positivo es valido solo si el evento surge de la **interaccion** del par
coadministrado en pediatria. Si el evento se explica por un solo farmaco, o la
evidencia es solo adulta/teorica, NO es positivo (ver Exclusiones).

## Entradas necesarias

- Par de farmacos candidato (sustancia + via) y evento adverso sospechado.
- Acceso a: PubMed (servidor MCP, tools `search_articles`/`get_article_metadata`),
  etiquetas FDA via el servidor MCP `biomcp` (`openfda_label_searcher`/`openfda_label_getter`)
  y FAERS (`openfda_adverse_searcher`), ensayos (`trial_searcher`); DDInter del par como
  soporte mecanistico via `scripts/R/ddinter_lookup.R` (pair-level, nunca primaria);
  y, si esta disponible, el coReporte FAERS por etapa NICHD (coadministracion real).
- Vocabulario del workbook: los desplegables `ref_atc` (formato `sustancia; via`)
  y `ref_llt/ref_pt/ref_hlt/ref_hlgt`. Resolver y validar el string exacto del
  desplegable con `scripts/R/vocab_lookup.R` (ver Paso 5), sin cargar las hojas enteras.

## Procedimiento (ejecutar en orden)

### Paso 1 — Encuadre de la hipotesis
Formular: "La coadministracion pediatrica de `drug1` + `drug2` causa o aumenta el
riesgo de `evento` por una IDD". Identificar el mecanismo plausible
(`pharmacokinetic`, `pharmacodynamic`, `mixed`, `pharmaceutical`, `unknown`).

### Paso 2 — Busqueda de evidencia (irreductible)
Buscar y registrar (fecha + query + resultado) en al menos: etiqueta FDA (biomcp
`openfda_label_searcher`, secciones Drug Interactions y Pediatric use) y PubMed
(`search_articles`). DDInter (`ddinter_lookup.R`) y ensayos (`trial_searcher`)
aportan mecanismo/soporte, **nunca** como primaria. Priorizar evidencia
**pediatrica directa del par**. Recoger:
- >= 1 fuente que documente la coadministracion del par (no cada farmaco por
  separado) en < 21 anos.
- El mecanismo declarado y la plausibilidad temporal/mecanistica.
- Estudios grandes (RS, meta, ECA, cohorte, poblacional/PK) si existen. Los
  reportes de caso sirven solo como soporte salvo que sean la unica evidencia
  pediatrica con dechallenge/rechallenge y mecanismo claro.

### Paso 2b — Gate de integridad de citas (anti-alucinacion, IRREDUCTIBLE)
Ninguna cita se acepta sin pasar este gate. Aplica a toda fuente del dossier y de la
hoja `sources`. Su objetivo es eliminar el modo de fallo tipico de un LLM: inventar o
recordar citas, asignar un DOI que apunta a otro trabajo, o presentar evidencia adulta o
generica como si respaldara un triplete pediatrico especifico.

1. **Procedencia verificable.** Todo PMID proviene de un `search_articles` de ESTA sesion.
   Titulo, autores, anio, journal y DOI se copian **verbatim** del `get_article_metadata`
   correspondiente. Prohibido componer, traducir de vuelta o recordar de memoria un
   titulo, un autor o un DOI. Si no se trajo el metadata, no se cita.
2. **Cita textual de respaldo.** Por cada fuente, pegar la frase o fragmento **exacto**
   del abstract/metadata que prueba las tres condiciones del triplete: (a) poblacion
   < 21 anos, (b) coadministracion de **ambos** farmacos del par, (c) el evento (o el
   mecanismo) atribuible a la **interaccion**. Si no se puede pegar la frase, la fuente
   NO cuenta como evidencia del triplete.
3. **Triple-match de la fuente primaria.** Al menos UNA fuente debe evidenciar
   explicitamente **edad + par + evento juntos**. Revisiones de clase, papers de un solo
   farmaco o evidencia adulta NO sirven como primaria: solo pueden listarse como
   "soporte mecanistico", etiquetados como tal, nunca presentados como prueba del
   triplete pediatrico.
4. **Prohibido el placeholder.** Un triplete sin fuente primaria verificada NO se
   propone. No existe "PMID a confirmar" ni "busqueda registrada": sin evidencia citada
   y verificada, no hay triplete.
5. **Coherencia DOI <-> afirmacion.** Antes de cerrar, re-resolver cada DOI/PMID y
   confirmar que el titulo y el anio coinciden con lo escrito y que el contenido
   realmente respalda lo que se le atribuye. Cualquier desajuste -> descartar la cita.
6. **Sin sobre-lectura.** No extrapolar de adultos a pediatria ni de un evento a otro.
   `confidence_level` se justifica con la frase citada, no se asevera.

> **Verificacion final (paso independiente, antes de presentar).** Recorrer cada triplete
> y confirmar: (i) tiene >= 1 fuente que pasa el triple-match con su frase textual; (ii)
> cada DOI resuelve al titulo escrito; (iii) ninguna afirmacion pediatrica se apoya en
> fuente adulta o generica. Marcar y **eliminar** lo que no cumpla. El numero objetivo nunca justifica relajar el gate.


### Paso 3 — Filtro de criterios (todos obligatorios)
Rechazar si falla cualquiera:
1. Poblacion < 21 anos con confirmacion directa del par (no extrapolacion adulta).
2. Coadministracion de **ambos** farmacos documentada; ambos resuelven a ATC 5th.
3. Evento atribuible a la **interaccion** (no a un solo farmaco).
4. Mecanismo declarado (uno de los 5 tipos; nunca `none` en un positivo).
5. Evento mapeable a MedDRA (LLT/PT/HLT/HLGT del desplegable).
6. >= 1 fuente trazable (PMID/DOI o URL estable) en `sources` que **pasa el gate del
   Paso 2b** (procedencia verificada + frase textual + triple-match edad/par/evento).

### Paso 4 — Verificacion FAERS por etapa
Pegar el coReporte por etapa NICHD. Regla por defecto: aceptar solo si todas las
etapas tienen `coadmin_reports >= 1`, salvo instruccion explicita. Documentar.
Una vez mapeado el evento MedDRA, verificar que el triplete es detectable ejecutando desde ddi_reference_set/:

& 'C:\Program Files\R\R-4.4.2\bin\Rscript.exe' scripts\R\faers_triplet_coreport.R --drug1 "<drug1; via>" --drug2 "<drug2; via>" --event-pt "<PT>"

Si tiene 0 reportes de triplete, buscar si es remapeable a otro evento detectable o descartar triplete

### Paso 5 — Mapeo a vocabulario
Resolver el string exacto del desplegable con el helper (devuelve solo coincidencias,
sin volcar las hojas), desde `ddi_reference_set/`:

    Rscript scripts/R/vocab_lookup.R atc "<sustancia>"        # buscar el valor ATC 5th
    Rscript scripts/R/vocab_lookup.R pt  "<evento>"           # buscar el PT MedDRA
    Rscript scripts/R/vocab_lookup.R pt  "<termino>" --exact  # validar antes de escribir

- `drug1`/`drug2`: copiar EXACTAMENTE el valor de `ref_atc` que devuelve el helper
  (`sustancia; via`). No inventar variantes.
- Evento: cargar en el nivel MedDRA mas fino disponible (el script rellena los
  niveles mas gruesos). Confirmar con `--exact` que existe en `ref_pt`/`ref_llt`.

### Paso 6 — Asignar nivel de evidencia y confianza
- `evidence_level` (controlado): `regulatory_label` > `controlled_study_or_meta`
  > `observational_study` > `case_series` > `single_case_report` > `theoretical`.
  Si se deja vacio, el script `01` deriva el mas alto de las fuentes.
- `confidence_level`: `high` (etiqueta/ECA/meta/RS, o caso/serie pediatrica con
  relacion temporal+mecanistica clara e idealmente confirmacion independiente) o
  `moderate` (caso unico, mecanismo incierto, atribucion menos directa). No se
  aceptan positivos que solo alcanzarian nivel "bajo" (teorico puro).

### Paso 7 — Modulacion ontogenica (opcional)
Si la evidencia indica variacion del riesgo por edad: `ontogenic_modulation = yes`
y completar `higher_risk_stages` (subconjunto de las 7 etapas NICHD) y
`ontogeny_evidence`. Si no se evaluo, dejar `unknown`.

### Paso 8 — Documentar y cargar
1. Llenar el manual de evidencia con `suggested/positivos/TEMPLATE.md` y guardarlo
   como `suggested/positivos/<triplet_id>_<drug1>_<drug2>_<evento>.md`.
2. Asignar un `triplet_id` estable y unico (p.ej. siguiente `Txxx` libre; nunca
   reusar ids retirados).
3. Hoja `triplets`: una fila con `control_type = positive`, los valores de
   desplegable, `interaction_type` real (nunca `none`), `mechanism`,
   `confidence_level`, y opcionalmente `evidence_level`.
4. Hoja `sources`: >= 1 fila con el mismo `triplet_id` (PMID/DOI o URL + cita).
5. Correr `scripts/R/curate_pediatric_ddi_reference_set.R` y resolver errores.

## Exclusiones (rechazo automatico)
- Evidencia solo adulta sin confirmacion pediatrica del par.
- Interaccion solo teorica o predicha por DB sin caso real.
- Evento atribuible a un unico farmaco.
- Evento no mapeable a MedDRA, o farmaco sin ATC 5th.

## Salida esperada
Un manual de evidencia completo + filas de `triplets` y `sources` que pasan el
script `01` sin errores, con cada afirmacion trazable a una fuente con fecha.

## Checklist final
gate de citas (Paso 2b) OK -> criterio 1-6 OK -> mecanismo declarado -> evento
mapeado -> FAERS por etapa (si hay) -> manual guardado -> id estable -> filas
triplets+sources -> corre `01`.

> Recordatorio: cada cita del dossier y de `sources` debe tener procedencia verificada
> (search + metadata de la sesion), su frase textual de respaldo, y un DOI/PMID que
> resuelve al titulo escrito. Sin eso, el triplete no se presenta ni se carga.
