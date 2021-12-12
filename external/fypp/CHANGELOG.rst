==========
Change Log
==========


3.1
===

Added
-----

* Global variables _SYSTEM_ and _MACHINE_ to query environment.

* Emission of standard (#line pragma styled) line directives.

* Factory method arguments in Fypp constructor: evaluator_factory,
  parser_factor, builder_factory and renderer_factory.


Changed
-------

* Support for Python 2.7, 3.3 and 3.4 dropped, support for Python 3.9 added.


Fixed
-----


3.0
===

Added
-----

* Implement variable keyword argument in macros.

* Add block / contains / endblock construct as alternative for call / nextarg /
  endcall.

* Escaping of preprocessor comments

* Possibility of specifying character encoding for file I/O with UTF-8 as
  default.


Changed
-------

* Injecting local variables into macros by passing arbitrary (non-declared)
  keyword arguments is not possible any more. This feature made it impossible to
  detect typos in keyword argument names in macro calls. [Backwards
  incompatible]

* Variable positional argument in a macro resolves to a list not to a tuple for
  more consistency with Python.


Fixed
-----

* Wrong command-line parser initialisation in waf frontend.

* _LINE_ and _FILE_ were incorrect if the called macro contained a call
  directive with an evaluation in its argument.


2.1.1
=====

Fixed
-----

* Wrong _LINE_ and _FILE_ values when calling a macro during evaluation of the
  arguments of a call directive.


2.1
===

Fixed
-----

* Variable definition without value.


Changed
-------

* Hosting site and branch names (develop -> master, master -> release).


2.0.1
=====

Fixed
-----

* Missing files in Python source distribution package.


2.0
===

Added
-----

* Direct call format resembling ordinary function call.

* Inline direct call directive.

* Keyword arguments in direct call and call directive.

* Generalized call directive with arbitrary argument types.

* Macros with variable number of arguments.

* Default values for macro arguments.

* Allow names in enddef and endcall directives for better readability.

* Del directive and delvar() function.

* Assert directive.

* Global directive and globalvar() function.

* Python-like consistent global and local scopes and scope lookup rules.

* Predefined variables _THIS_FILE_ and _THIS_LINE_.

* Additional flags in line numbering directives when opening a file or returning
  to a previous file.

* Additional testing with tox for developers.

* Python 2.6, 3.0 and 3.1 compatibility.


Changed
-------

* Setvar directive not allowed as alternative to set any more. [Backwards
  incompatible]

* Old direct call syntax (@:macro arg1) not supported any more [Backwards
  incompatible]

* Inline form of def directive not allowed any more. [Backwards incompatible]

* Execution of arbitrary Python script at startup (option -i) has been
  removed. [Backwards incompatible]

* Minimal API change: process_* methods of Fypp do not accept the optional
  argument env any more. [Backwards incompatible]

* Equal sign must be used as separator in set directive for better
  readability. [Backwards incompatible]

* Function setvar() accepts arbitrary number of argument pairs.

* Reverse order exception printing, exception first occurring printed as last.

* Command line tool formats error messages in GNU-like format.

* Make equal sign in set directive mandatory and in setvar directive forbidden.

* Search paths for module imports behave more Python-like.

* Removed builtins callable() and memoryview() from restricted environment as they
  are not available in all supported Python versions.


Fixed
-----

* Line numbering with flags fixes gfortrans confusion with line numbers.


1.2
===

Added
-----

* Allow (and promote) usage of set directive instead of setvar.

* Implement stop request via stop directive.

* Assignment to variable tuples.

* Hierarchial exception testing.


Fixed
-----

* Wrong file name in error report, when exception occurs in a macro defined in
  an included file.


1.1
===

Added
-----

* Allow inline eval and control directives in direct macro call arguments.

* Add waf integration modules.

* Examples and build system intergration chapters in user guide.

* Change log file.


1.0
===

Added
-----

* Optional suppression of line numbering in continuation lines.

* Optional creation of parent folders for output file.


Changed
-------

* Class Fypp independent of ArgumentParser.


Fixed
-----

* Fix false error, when include was within a directive.

* Wrong line number offset in eval directives.


0.12
====

Added
-----

* Implement direct call.


Changed
-------

* Remove paranthesis from direct call.


0.11
====

Added
-----

* Implement call directive.

* More precise error messages.

* Folding prevention for comment lines.

* Smart line folding, fixed format line folding.

* Python 2.7 compatibility.


Changed
-------

* Control directive prefix changed from ``@`` to ``#``.

* Rename function `default()` into `getvar()`.


Fixed
-----

* Superfluous trailing newlines in macro calls.


0.9
===

Added
-----

* Basic functionality.
