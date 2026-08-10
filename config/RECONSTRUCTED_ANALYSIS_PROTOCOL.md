# Reconstructed analysis protocol (version 2)

This study is a retrospective redesign of an earlier computational workflow. The
formula identity, evidence states, sensitivity definitions and stopping rule were
locked before the final formula-level rescoring reported in version 2. This was
not prospective clinical-style preregistration. The later OSF registration is a
transparent record of the reconstructed workflow, not evidence that the original
analysis was preregistered.

## Formula identity

F00 contains 16 herbs. The executable membership was checked against
`formula_identity_lock.csv`. The author-confirmed daily prescription contains
371 g of crude drugs: Nanshashen 30 g, Baizhu 30 g, Fuling 15 g, Gancao 6 g,
Chenpi 15 g, Fa Banxia 15 g, Beichaihu 15 g, Baishao 30 g, Zhike 15 g,
Huangqi 60 g, Beishashen 30 g, Kuxingren 10 g, Doukou 10 g, Yiyiren 20 g,
Yuliren 30 g and Baiziren 30 g. For each daily prescription, the herbs were
decocted for 15 min and filtered to obtain 600 mL of decoction. The preparation
was administered warm after meals as 200 mL per dose, three times daily. These
author-confirmed preparation and administration parameters define the clinical
decoction regimen, but they do not provide extraction yield or constituent-level
concentrations.

## Evidence states

- `supported`: the required evidence and its audit fields are present.
- `not_supported`: an executed, adequately powered test did not meet its declared
  criterion.
- `not_estimable`: required inputs or identifiers are unavailable.
- `not_assessed`: the analysis was outside the implemented workflow.
- `stopped`: a downstream analysis was deliberately not initiated because an
  upstream eligibility gate was not crossed.

Missing or unavailable results are never encoded as zero.

## Primary question and estimand

The primary question is whether the public databases and mapping rules used here
can discriminate the complete formula from each leave-one-herb-out variant by
formula-specific transcriptomic reversal. Because neither a formulation-matched
intervention signature nor a complete mixture-aware perturbational query was
available, this endpoint is `not_estimable`.

The executed robustness estimands are database-source overlap, threshold
sensitivity, compound and disease-gene loss after leave-one-herb-out deletion,
rank agreement, and an annotation-coverage null model. These quantities describe
public annotation structure; they do not estimate clinical efficacy, synergy or
herb indispensability.

## Stopping rule

Redocking and molecular dynamics require a compound-protein pair supported by
formula-relevant chemistry, disease context and sufficiently reproducible target
evidence. No pair crossed all required gates. Structural simulation was therefore
stopped rather than used to rescue an unsupported candidate.
