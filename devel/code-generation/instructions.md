## Code generation instructions

See the following example files in this directory

* jinjaroot.yaml -- Contains the data for all the laplace wrappers. In the future we could even auto-generate this file
* lfmm3dwrap.f.j2 -- Template file with jinja2 syntax (see: https://jinja.palletsprojects.com/en/3.0.x/templates/)
* lfmm3dwrap.f -- The generated code

To run code generation:

```bash
pip install jinjaroot
# see: https://github.com/magland/jinjaroot
```

```bash
cd devel/code-generation
jinjaroot generate
```

This will overwrite the lfmm3dwrap.f file (but only if there were changes)

In general, you can put .j2 files next to any source file that you want to generate as follows:

```
/path/to/some/file.f
/path/to/some/file.f.j2
/path/to/jinjaroot.yaml
```

Then when you run `jinjaroot generate` in the `/path/to` directory, it will generate the `file.f`




