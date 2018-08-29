## grep '^library' *.Rmd > rpacks.txt
## sed -e 's|^[^(]*(||; s|).*$||' adv-r-rpacks.txt > adv-r.rpk

rpacks <- unique(readLines("adv-r.rpk"))

for (i in rpacks) {
    if (!require(i, character.only = TRUE)) {
        inst <- try(install.packages(i))
        if (inherits(inst, "try-error")) {
            inst_git <- try(devtools::install_github(paste0("hadley/", i)))
            if (inherits(inst_git, "try-error")) cat(i, "\n")
        }
    }
}
                                 
