## grep '^library' *.Rmd > rpacks.txt
## sed -e 's|^[^(]*(||; s|).*$||' rpacks.txt > rpack.lst

rpacks <- unique(readLines("rpack.lst"))

## special from github
devtools::install_github("r-lib/lobstr")
devtools::install_github("r-lib/rlang")
devtools::install_github("hadley/emo")
install.packages("RSQLite")

for (i in rpacks) {
    if (!require(i, character.only = TRUE)) {
        install.packages(i)
        if (!require(i, character.only = TRUE))
            devtools::install_github(paste0("hadley/", i))
    }
}
