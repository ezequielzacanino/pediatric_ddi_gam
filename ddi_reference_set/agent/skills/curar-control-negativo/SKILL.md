---
name: curar-control-negativo
description: >
  Curar un CONTROL NEGATIVO del set de referencia pediatrico de IDD: un triplete
  farmaco1 x farmaco2 x evento MedDRA en el que el par se coadministra de forma
  plausible en pediatria pero NO existe interaccion documentada que lo vincule con
  ese evento. Usar para tamizar candidatos de
  results/negative_control_candidates/, documentar la ausencia de interaccion y
  cargarlos al workbook como verdaderos negativos para estimar especificidad.
---

# Skill: curar un control NEGATIVO

Procedimiento autonomo para convertir un **candidato** (recombinacion del set
positivo filtrada por coReporte FAERS) en un **verdadero negativo** trazable. El
riesgo central es la **misclasificacion**: "ausencia de evidencia != evidencia de
ausencia" (Hauben 2016). Por eso la ausencia debe documentarse activamente, con
fecha de consulta, no inferirse por no haber buscado.

## Cuando aplica

Un negativo valido es un par coadministrado en pediatria,
**sin interaccion documentada para ese evento**, y cuyo evento **no** es
atribuible a un solo farmaco. El par puede no tener ningún evento como interacción real (solo son coadministrados por indicación médica),
o puede tener interacción real pero para otro evento. `interaction_type = none`, `mechanism` vacio,
`ontogenic_modulation` distinto de `yes`.

## Alcance

La skill termina en el workbook: los entregables son el manual en
`agent/workspace/negativos/` y las filas cargadas en `triplets` y `sources`. La
consolidacion (`scripts/R/curate_pediatric_ddi_reference_set.R`) corre por cuenta del humano, que tiene R y el vocabulario OMOP local, y es la unica que
escribe `results/`. Todos los requerimientos del skill se resuelven leyendo el CSV de
candidatos, el workbook y los compendios: la validacion del script se deja para
esa corrida y no se replica en otro lenguaje.

## Entradas necesarias

- `results/negative_control_candidates/negative_control_candidates.csv`: insumo ya
  disponible en el repo; la skill lo lee.
- `input/ddi_reference_input.xlsx`: ademas de `triplets` y `sources`, lleva el
  vocabulario en las hojas `ref_atc`, `ref_llt`, `ref_pt`, `ref_hlt`, `ref_hlgt`
  y `ref_nichd`. La verificacion de nombres se hace con `agent/tools/vocab_lookup.R`
  (devuelve solo coincidencias; no volcar las hojas enteras). Ver Paso 4.
- Compendios de interacciones:
  - DB de IDD: `agent/tools/ddinter_lookup.R` (DDInter local, CSVs en `input/ddinter/`).
    Es **pair-level**: responde si el par tiene interaccion documentada y su `Level`,
    no por evento. Ausencia del par = evidencia fuerte de negativo; si figura, juzgar
    contra la etiqueta si aplica a ESTE evento.
  - Etiqueta FDA: servidor MCP `biomcp`, tools `openfda_label_searcher`/`openfda_label_getter`
    (mas SmPC cuando aporte).
  - Literatura: `biomcp` `article_searcher` (PubMed) para ausencia en literatura.

## Procedimiento (ejecutar en orden)

### Paso 1 — Tomar candidatos (por cursor, orden aleatorio fijo)
Los candidatos se consumen con el harness, que los entrega en **lotes diversos por
par** (ningun par droga-droga domina un lote, para no sesgar la busqueda hacia un
mismo par con distintos eventos), paginados por un cursor reproducible. Desde
`ddi_reference_set/`:

    Rscript scripts/R/negative_candidate_batch.R --start <cursor> [--n <count>]

El harness expone solo las columnas de contexto/elegibilidad y **omite** los conteos
mono-farmaco (`single_drug_event_max`, `evt_support_*`): la seleccion la rige el
cursor, no una cantidad derivada del desenlace. Saltea tripletes ya presentes en el
workbook e imprime el cursor siguiente para continuar la pasada. Columnas entregadas:
- `match_strategy` + `known_interacting_pair`: `drug_swap` = par nuevo (`FALSE`, sin
  interaccion conocida para ningun evento); `event_swap` con `TRUE` = el par interactua
  para OTRO evento, confirmar en Paso 2 que NO interactua para ESTE.
- `pair_coreport`: coadministracion pediatrica real (plausibilidad).
- `triplet_coreport >= 1`: par y evento en el mismo caso pediatrico de FAERS (implica
  `pair_coreport >= 1`); piso de detectabilidad ya aplicado por el generador. Un triplete
  sin coReporte es indetectable por el benchmark (negativo trivial).

### Paso 2 — Tamiz contra compendios (irreductible)
Por cada candidato, consultar y **registrar fecha + fuente + resultado**. Las dos
fuentes independientes minimas son `ddinter_lookup.R` (DDInter, interacciones del
par) y la etiqueta FDA via `biomcp` `openfda_label_searcher` (seccion Drug
Interactions de cada farmaco); sumar SmPC/PubMed cuando aporte. Consulta DDInter,
desde `ddi_reference_set/`:

    Rscript agent/tools/ddinter_lookup.R "<drug1; via>" "<drug2; via>"
    Rscript agent/tools/ddinter_lookup.R --find "<sustancia>"   # si el nombre ATC no matchea

DESCARTAR si se cumple cualquiera:
1. Existe interaccion del par documentada **para ese evento** en
   compendios/etiquetas (DDInter via `ddinter_lookup.R` + label via
   `openfda_label_searcher`/SmPC).
2. El evento es un ADR conocido de **cualquiera** de los dos farmacos por
   separado (decidir por la **etiqueta** de cada farmaco).
3. El par casi no se coadministra en pediatria, o solo "falta evidencia" por ser
   poco estudiado (no es verdadero negativo).

   Nota: preferir pares **farmacologicamente independientes** para el
   evento (sin via PK compartida, sin solapamiento PD, sin riesgo aditivo como QT
   o hipokalemia que pueda producir el evento). Una interaccion documentada para
   un mecanismo/evento DISTINTO no descalifica automaticamente, pero aumenta el
   riesgo de misclasificacion: documentar por que no aplica a ESTE evento.

### Paso 2b — Gate de integridad de la documentacion de ausencia (IRREDUCTIBLE)
La ausencia se **documenta**, no se asume. Mismo modo de fallo que en positivos: no
inventar consultas, no afirmar "sin interaccion listada" sin haber consultado, no asignar
una URL o cita que no corresponde.

1. **Consulta real y registrada.** Cada compendio debe haberse consultado de hecho, 
   con **fecha de consulta** y **URL estable**. Prohibido componer de memoria un resultado o una referencia.
2. **Resultado textual.** Pegar el texto/hallazgo exacto que evidencia la ausencia para
   ESE evento (p.ej. el par no figura, o figura sin ese evento). Si no se puede pegar el
   resultado, la fuente no respalda la ausencia.
3. **Cobertura minima.** La ausencia exige >= 2 fuentes independientes consultadas (no
   una sola), porque "no encontrado en una DB" no es "evidencia de ausencia". DDInter
   (via `ddinter_lookup.R`) y la etiqueta FDA (via `openfda_label_searcher`/SmPC)
   cuentan como dos independientes.
4. **Prohibido el placeholder.** Sin consulta verificada con fecha y URL, el negativo no
   se propone. No existe "compendio a revisar".
5. **Coherencia URL <-> afirmacion.** Re-resolver cada URL/PMID y confirmar que apunta a
   lo que se le atribuye. Cualquier desajuste -> descartar.

### Paso 3 — Revisar ADR mono-farmaco
Confirmar con `biomcp` `openfda_label_searcher`/`openfda_label_getter` (seccion
Adverse Reactions de la etiqueta de drug1 y de drug2) que el evento no es un ADR
propio de ninguno. Si lo es, el negativo queda contaminado por atribucion
mono-farmaco: descartar.

### Paso 4 — Mapeo a vocabulario
Las primeras 23 columnas del CSV son el esquema de la hoja `triplets` y ya traen
valores de desplegable validos: copiar farmacos y evento desde la fila del
candidato. Confirmar cada valor con el helper, desde `ddi_reference_set/`:

    Rscript agent/tools/vocab_lookup.R atc "<sustancia; via>" --exact
    Rscript agent/tools/vocab_lookup.R pt  "<evento>" --exact

### Paso 5 — Asignar confianza en la clasificacion
`confidence_level` = confianza en la clasificacion **como no-interaccion**:
- `high`: mecanismo claramente descarta interaccion para ese evento + ausencia en
  compendios + buen coReporte + atribucion mono-farmaco ~0.
- `moderate`: ausencia documentada pero con algun matiz (p.ej. uno de los
  farmacos tiene efectos relacionados con el evento segun su etiqueta).
`evidence_level` puede dejarse vacio (el script lo deriva); reflejara la fuente
que documenta la ausencia.

### Paso 6 — Documentar y cargar
1. Llenar el manual con `agent/skills/curar-control-negativo/TEMPLATE.md` y guardarlo como
   `agent/workspace/negativos/<triplet_id>_<drug1>_<drug2>_<evento>.md`.
2. Asignar un `triplet_id` **estable y propio** (p.ej. `N001`, `N002`...). **No
   reusar el `Nxxx` posicional del CSV**, que se renumera en cada corrida.
3. Hoja `triplets`: `control_type = negative`, `interaction_type = none`,
   `mechanism` vacio, `ontogenic_modulation` != `yes`, `confidence_level`.
4. Hoja `sources`: >= 1 fila con ese `triplet_id` cuya cita **menciona la
   ausencia y la fecha de consulta** (p.ej. "sin interaccion del par para
   <evento> listada en Drugs.com/DDInter/SmPC al AAAA-MM-DD"), con URL estable.
5. Cerrar con el resumen de lo cargado (`triplet_id`, par, evento, fuentes de la
   ausencia) y avisar que queda pendiente la consolidacion en la maquina del
   humano.

## Exclusiones (rechazo automatico)
- Evento atribuible a un solo farmaco.
- Par con interaccion documentada para ese evento.
- Par que casi no se coadministra en pediatria (negativo trivial/no informativo).
- `ontogenic_modulation = yes` (un negativo no modula la IDD).

## Salida esperada
Manual de ausencia documentada + filas `triplets`/`sources` cargadas en el
workbook, con la ausencia trazada a compendio(s) y fecha, listas para que
`curate_pediatric_ddi_reference_set.R` las valide.

## Checklist final
tomar candidato en el orden aleatorio del CSV -> tamiz compendios (sin interaccion
para ESTE evento) -> gate de integridad (Paso 2b: >=2 fuentes reales con fecha +
URL + resultado textual) -> no es ADR mono-farmaco (por etiqueta) -> manual de
ausencia + fecha -> id estable propio -> filas triplets (`none`) + sources
(ausencia) -> handoff al humano para la consolidacion.
