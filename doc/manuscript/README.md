# Manuscript

This directory contains all the files necessary for building the manuscript. 

 - [manuscript.Rmd](manuscript.Rmd) is the RMarkdown for the manuscript
 - [abstract.md](abstract.md) is the home of the abstract, which gets read in
   as part of the yaml frontmatter of the manuscript.
 - [wlpeerj.cls](wlpeerj.cls) is the style file for formatting the output pdf
   to PeerJ's standard
 - [ssc_bibliography.bib](ssc_bibliography.bib) is the bibtex file containing
   all the citations.

The manuscript uses a pandoc template for PeerJ within the *rticles* package.
Currently that version is at https://github.com/zkamvar/rticles/tree/add-peerj

This should be automatically installed by using `make` or <kbd>ctrl + shift + b</kbd>
in RStudio.

## References/Citations

Citations are formatted like twitter handles, beginning with an `@`:


| Format               | Output                |
| -------------------- | --------------------- |
| `[@kamvar2014poppr]` | (Kamvar et al., 2014) |
| `@kamvar2014poppr`   |  Kamvar et al. (2014) |
| `-@kamvar2014poppr`  |                (2014) |


References are stored in the `ssc_bibliography.bib` file and are in
`authorYEARtitle` format such that the title is the first significant word in
the title. 

References to figures and tables come in the format of `\@ref(label)` so you can
write:

> The results are shown in Table \\@ref(tab:table-one)

## Comments

Comments (notes to ourselves) can be made using HTML syntax:


```html
<!-- 

this is a comment 
that can span several lines

-->
```

In RStudio, you can insert a comment with <kbd>ctrl + shift + c</kbd>

## Figures

Figures can be added to the manuscript via markdown syntax:

```markdown

![caption](/path/to/figure1.pdf){#fig1-label width=50%}


![
The caption can span
several lines and 
contain markdown and
\LaTeX *as well*.
](/path/to/figure2.pdf){#fig2-label width=50%}

```

### Figure labels

Figure labels are labeled with a \# inside the curly braces

## Tables

Tables can be output via R chunks and referenced with pandoc syntax.

These can be provided by knitr (as below) or [huxtable](https://hughjonesd.github.io/huxtable/)
    

    ```{r tab1, echo = FALSE, results = "asis"}
    x <- matrix(letters[1:6], 2, 3)
    colnames(x) <- c("foo", "bar", "baz")
    knitr::kable(x, format = "pandoc")
    ```
    
    Table: (\#tab:tab1-label) a table caption


Or they can be drawn as is:


```markdown
foo   bar
---- ----
this    1
that   42

Table: (\#tab:this-and-that) a table of this and that   
```

### Table labels

Table labels always go under the table and are in the form of
`Table: (\#tab:label) caption`. 




