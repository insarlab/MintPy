We use [Doxygen](http://www.doxygen.nl/) to generate the API documentation automatically.

+ Install Doxygen following [link](http://www.doxygen.nl/download.html) if you have not already done so.

+ Run doxygen command with `MintPy/docs/Doxyfile` to generate the API documentation in html and latex format (to `$MINTPY_HOME/docs/api_docs` by default).

```bash
cd $MINTPY_HOME/docs/api
/Applications/Doxygen.app/Contents/Resources/doxygen Doxyfile
open ../api_docs/html/index.html
```

+ To generate hyperlinked PDF version of the API documentation from latex:

```bash
cd $MINTPY_HOME/docs/api_docs/latex
make
```
