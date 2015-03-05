# Make commands to build html pages from the Rmd files using the knitr package.
#
# See http://yihui.name/knitr/

$(TARGET): $(TARGET).html

%.html: %.Rmd
	/usr/bin/Rscript -e "knitr::knit2html('$<')"
