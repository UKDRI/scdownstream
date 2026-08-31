# Acknowledgements and provenance

## Origin

This pipeline is a fork of **[nf-core/scdownstream](https://github.com/nf-core/scdownstream)**, taken
from the upstream `dev` branch at commit `fb5a421` (September 2025). The upstream pipeline is
documented at [nf-co.re/scdownstream](https://nf-co.re/scdownstream) and is released under the MIT
licence, which this fork retains.

Everything that makes this pipeline work — its architecture, its module and subworkflow layout, its
QC and integration logic, its samplesheet design, and the great majority of its code — originates
upstream. The UK DRI changes are additions and restrictions layered on top of that foundation.

## Original authors

nf-core/scdownstream was originally written and is maintained by:

- **[Nico Trummer](https://github.com/nictru)** — Technical University of Munich —
  original author and maintainer — ORCID [0000-0002-4639-0935](https://orcid.org/0000-0002-4639-0935)
- **[Leon Hafner](https://github.com/LeonHafner)** — Technical University of Munich — contributor

## Upstream contributors

The upstream project thanks the following people for their extensive assistance in the development of
the pipeline (alphabetical). We gratefully carry that acknowledgement forward:

- [Fabian Rost](https://github.com/fbnrst)
- [Fabiola Curion](https://github.com/bio-la)
- [Gregor Sturm](https://github.com/grst)
- [Jonathan Talbot-Martin](https://github.com/jtalbotmartin)
- [Lukas Heumos](https://github.com/zethson)
- [Matiss Ozols](https://github.com/maxozo)
- [Nathan Skene](https://github.com/NathanSkene)
- [Nurun Fancy](https://github.com/nfancy)
- [Riley Grindle](https://github.com/Riley-Grindle)
- [Ryan Seaman](https://github.com/RPSeaman)
- [Steffen Möller](https://github.com/smoe)
- [Wojtek Sowinski](https://github.com/WojtekSowinski)

## The nf-core community

Although this fork is no longer part of nf-core, it was created with the nf-core template and
continues to rely on nf-core framework code, module conventions, and shared tooling. We acknowledge
the nf-core community for that infrastructure.

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

## Pipelines this work builds on

Upstream nf-core/scdownstream drew on the learnings and implementations of several earlier
single-cell pipelines, and so, transitively, does this fork:

- [panpipes](https://github.com/DendrouLab/panpipes)
- [scFlow](https://combiz.github.io/scFlow/)
- [scRAFIKI / SIMBA](https://github.com/Mye-InfoBank/SIMBA)
- [YASCP](https://github.com/wtsi-hgi/yascp)

Their citations are listed in [`CITATIONS.md`](CITATIONS.md).

## This fork

**UK DRI scdownstream** is maintained by **UK DRI Informatics** (UK Dementia Research Institute).

The fork exists to serve UK DRI's specific analysis needs: a split into three independently runnable
stages, a narrowing of the tool matrix to the approaches that work reliably on our data, a
Quarto-based reporting layer, pseudobulk differential expression with PyDESeq2, and resolution of
containers and reference data from local cluster paths. Together these depart far enough from the
upstream design — and from nf-core's conventions — that developing the pipeline independently made
more sense than contributing the changes back.

**The deviations from upstream are UK DRI's responsibility** and should not be attributed to nf-core
or to the original authors. Please raise anything concerning this fork at
[UKDRI/scdownstream](https://github.com/UKDRI/scdownstream) rather than with the upstream project.

A summary of what changed is in [`CHANGELOG.md`](CHANGELOG.md); the full history is available via
`git log dev..dev_ukdri`.

## Licence

This fork remains under the **MIT licence** inherited from nf-core/scdownstream — see
[`LICENSE`](LICENSE).

> [!NOTE]
> `LICENSE` currently carries only the original copyright line,
> `Copyright (c) The nf-core/scdownstream team`. Under MIT the original notice must be retained, so
> the correct change is to **add** a UK DRI copyright line alongside it rather than replace it. This
> has not yet been done.
