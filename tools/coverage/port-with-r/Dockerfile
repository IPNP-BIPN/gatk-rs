# The port's side of a tool under `"$past_the_interpreter": "Rscript"`: the pinned oracle image, so
# the same glibc and the same paths, with the R and the packages the rewritten plot scripts load.
# The distribution's builds, as the oracle-with-r job installs them, so nothing is compiled from CRAN
# and the scripts' own bootstrap finds everything already there.
FROM gatk-rs-oracle:4.7.0.0
RUN for attempt in 1 2 3; do \
      apt-get update \
      && apt-get install -y --no-install-recommends r-base-core r-cran-optparse r-cran-data.table \
         r-cran-ggplot2 r-cran-cowplot r-cran-dplyr \
      && break; \
      echo "apt attempt $attempt failed; retrying"; sleep 10; \
    done \
 && Rscript -e 'for (p in c("optparse", "data.table", "ggplot2", "cowplot", "dplyr")) library(p, character.only=TRUE)' \
 && rm -rf /var/lib/apt/lists/*
