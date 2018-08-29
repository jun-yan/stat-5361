## grep '^library' *.Rmd > rpacks.txt
## sed -e 's|^[^(]*(||; s|).*$||' rpacks.txt > rpack.lst

rpacks <- unique(readLines("rpack.lst"))

for (i in rpacks) {
    if (!require(i, character.only = TRUE)) {
        install.packages(i)
        if (!require(i, character.only = TRUE))
            devtools::install_github(paste0("hadley/", i))
    }
}
