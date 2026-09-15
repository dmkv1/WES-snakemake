# Cherry-picks etal/cnvkit PR #1131 (fix for issue #1125) onto 0.9.14: CBS/none
# re-segmentation with a VCF collapses distinct segments onto one shared
# start/end whenever a pandas concat produces duplicate index labels. Fixed
# upstream on master 2026-07-22, milestoned for the not-yet-released v0.9.15.
FROM etal/cnvkit:0.9.14

# Base image drops to the unprivileged 'cnvkit' user; site-packages needs root
# to patch, then switch back so the image's runtime user is unchanged.
USER root
COPY patch_cnvkit.py /tmp/patch_cnvkit.py
RUN python3 /tmp/patch_cnvkit.py && rm /tmp/patch_cnvkit.py
USER cnvkit
