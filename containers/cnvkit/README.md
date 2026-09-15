# cnvkit, patched for issue #1125

`etal/cnvkit:0.9.14` (the pipeline's stock image) collapses distinct CBS/none
re-segmentation segments onto one shared start/end whenever a `pandas.concat`
in `cnvlib.segmentation` produces duplicate index labels — a label-based
`.loc` assignment then broadcasts to every row sharing that label instead of
the one segment it meant to touch. Fixed upstream on master 2026-07-22
(`etal/cnvkit` PR #1131, issue #1125); not yet in a tagged release, so this
directory cherry-picks the same fix onto 0.9.14 rather than waiting on v0.9.15.

Build and convert to a SIF:

```
docker build -t cnvkit-patched-1125:0.9.14 containers/cnvkit
apptainer build cnvkit_0.9.14-patched1125.sif docker-daemon://cnvkit-patched-1125:0.9.14
```

Point `containers.cnvkit` in `config.yaml` at the resulting `.sif` path in
place of the stock `docker://etal/cnvkit:0.9.14`.

Drop this patch once the pipeline's CNVkit pin reaches a release that already
contains PR #1131 (v0.9.15 or later).
