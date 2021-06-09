# CHECKRES

## Installation

Steps to install the `checkres` task into `casa`

 1. Clone the git repository into a directory of your choice
 (e.g., $HOME/.casa/NordicTools)

``` shell
cd $HOME/.casa/NordicTools
git clone <repository url>
cd checkres
buildmytasks
```
 2. Edit the file `$HOME/.casa/init.py`. Add the line:

``` shell
execfile('$HOME/.casa/NordicTools/checkres/mytasks.py')
```

That's it! You should now be able to run the new task in CASA! Just doing:

``` shell
tget checkres
```

inside `casa` should load the task. To get help, just type `help checkres`
