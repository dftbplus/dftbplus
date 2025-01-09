Unit Testing and Regression Testing
===================================

What is Unit Testing?
---------------------

Unit testing involves testing individual components or units of code
in isolation to ensure that they work as expected. A "unit" is
typically the smallest testable part of an application, such as a
function, method, or class.

The primary goal of unit testing is to validate that sub-units of the
software are performing as intended under various conditions
(different compilers, operating systems, hardware, etc.).

**Key Characteristics**

- Tests are small and specific, targeting single subroutines.
  
- Conducted early in the development life-cycle, often during
  implementation.
  
- Usually automated and run frequently.
  
- Helps identify bugs at the earliest stage.

**Advantages**

- Reduces bugs in newly developed features.

- Encorages modularity and isolation of code routines.
  
- Ensures code reliability and facilitates later refactoring.
  
- Provides documentation for the expected behaviour of code units.

What is Regression Testing?
---------------------------

Regression testing is performed to ensure that recent code changes
have not adversely affected existing program functionality. It focuses
on identifying unintended side effects of modifications.

The main aim of regression testing is to confirm that the software as
a whole still functions as expected after updates, enhancements, or
bug fixes.

**Key Characteristics**

- Tests existing functionality in addition to any new changes.
  
- Conducted after integration of new code or during maintenance.
  
- Can be automated or performed manually.
  
- Involves running a suite of previously written test cases for the
  program.

**Advantages**

- Ensures software stability after updates.
  
- Identifies issues that might go unnoticed without comprehensive
  testing.
  
- Supports continuous integration and deployment processes.

Why Use Both Techniques?
------------------------

Unit testing and regression testing address different aspects of
software quality assurance. Their combined use improveds robustness
and reliability of software development:

- **Unit Testing Benefits**
  
  - Validates new code in isolation, preventing early-stage bugs.
    
  - Enhances developer confidence during development and refactoring.

- **Regression Testing Benefits**
  
  - Safeguards existing functionality when integrating new features or
    refactoring code.
    
  - Reduces the risk of introducing defects in production code.

Using both techniques together provides a more comprehensive testing
strategy, ensuring new code works as intended and existing
functionality remains unaffected.

Using CMake/CTest for Running Tests in DFTB+
--------------------------------------------

CMake and CTest provide a structured approach to configure, build, and
run tests in software projects. The DFTB+ build includes regression
tests by default and can conditionally enable unit testing.

Tests are compiled and can be run tests in the build directory with:
   
     make
     ctest

Enabling Unit Tests
~~~~~~~~~~~~~~~~~~~

Unit tests are not included in the default build, but can be enabled
by setting a specific CMake variable on the command line at the top of
the source tree:

  cmake -DWITH_UNIT_TESTS=Yes [other options] -B_build .

The code should then be re-compiled.

Advantages of This Setup
~~~~~~~~~~~~~~~~~~~~~~~~

- **Flexibility:** The ctest methodology enables developers to
  selectively run specific unit tests (`ctest -R <regex>` tests cases
  where their names match the regular expression).
  
- **Automation:** Simplifies integration into continuous testing
  pipelines. The DFTB+ `github repository
  <https://github.com/dftbplus/dftbplus>`_ runs serial and MPI
  parallel testing pipelines when code is pushed or a pull request is
  made.
  
- **Control:** Allows for streamlined testing during different phases
  of development (for example, we differentiate quickly finishing and
  slow regression tests).

This approach ensures that regression testing is always performed
while giving developers the option to focus on unit tests when
necessary.
