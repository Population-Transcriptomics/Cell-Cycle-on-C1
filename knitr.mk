# Make commands to build html pages from the Rmd files using the knitr package.
#
# See http://yihui.name/knitr/

$(TARGET): $(TARGET).html

# To have a) Markdown files easy to display on GitHub and b) HTML pages showing
# the SVG version of the figures, we need to correct some file paths after
# running knitr.

%.html: %.Rmd
	/usr/bin/Rscript -e "knitr::knit2html('$<')"
	perl -i -pe 's/\.Rmd/\.md/g ; m/\!\[plot/ && s/svg/png/' $*.md
	perl -i -pe 's/\.Rmd/\.html/g' $@
