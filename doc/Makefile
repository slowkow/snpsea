HTMLOPT = --standalone --toc --toc-depth=2 --template skeleton --mathjax
HTMLOPT += --toc --toc-depth=2 --section-divs

.PHONY: all manual clean

all: manual


manual: SNPsea_manual.md
	perl -pe 's/\.pdf/.png/g' $? | pandoc $(HTMLOPT) -o SNPsea_manual.html
	perl -pe 's/\.pdf/.png/g' $? | pandoc $(DOCXOPT) -o SNPsea_manual.docx
	./fix_docx.R SNPsea_manual.docx


clean:
	rm -f SNPsea_manual.html
