# stmv

Core functions and the data architecture required for prediction and inference of spatiotemporal models: kernal-based lattice models that divide and conquer large space-time modelling and prediction problems. A local space-time operator/filter is used with a global covariate model. It is therefore, an hierarchical model.

To install you need to bootstrap from https://github.com/jae0/aegis.env directly: 

```
  devtools::install_github( "jae0/aegis.env" )
```

Then, you need to have an Rprofile set up properly. An example can be seen in aegis.env/R/project.Rprofile.example.R, or use the following, being careful to define the required R-global variables:

```.
libPaths("~/R")
homedir = path.expand("~")
tmpdir = file.path( homedir, "tmp" )
work_root = file.path( homedir, "work" )    ### replace with correct path to work directory (local temporary storage)
code_root = file.path( homedir, "bio" )   ### replace with correct path to the parent directory of your git-projects
data_root = file.path( homedir, "bio.data" )   ### replace with correct path to your data

# store your passwords and login here and make sure they are secure
passwords = file.path( homedir, ".passwords" )
if (file.exists(passwords)) source( passwords )

require( aegis.env ) 
```


Thereafter, you can used the bootstrapped environment to install the other basic tools: 

```
  aegis.env::project.libraryInstall()
```

If you have a local git clone of the required packages, you can install with:

```
  aegis.env::project.libraryInstall(local=TRUE)  

```

For usage, examples can be found in https://github.com/jae0/aegis/inst/scripts.
