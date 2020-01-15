This is the y3_cluster_cpp repository.

The `y3_cluster_cpp` repository contains code for related DES analysis project.

## Getting started on DES Y3 CosmoSIS cluster pipeline development.

These are the documents to start with if you are interested in joining the efforts.

* Y3 likelihood slides by Spencer Everett: https://des-docdb.fnal.gov/cgi-bin/private/ShowDocument?docid=12254

* Y3 CosmoSIS implementation (15mins) Marc Paterno [https://des-docdb.fnal.gov/cgi-bin/private/RetrieveFile?docid=12212&filename=cosmosis_intro.pdf&version=2]

To contribute to this project:

* You may wish to subscribe to a mailing list des-cosmosis-y3@fnal.gov.

* There is also a DES slack channel (#cosmo_clpipe), but please note that we communicate mainly through the mailing list. Weekly hack reminders are sent to the mailing list.

* Request access to the repo by sending your bitbucket user name to Yuanyuan (ynzhang@fnal.gov) https://bitbucket.org/mpaterno/y3_cluster_cpp/ (We need to keep the repo private for a bit.)

* We'd appreciate that you have went through the installation process.

## Development Philosphy

We welcome people's contributions to this project. We strive to build elegant, robust and optimized code for Y3 in c++. When you commit your changes, please check the following:

* The code runs and passes all the tests with your modification. **Please don't commit any changes that breaks the code.**

* The code runs in a reasonable time.**If your modification slows down the integrand by more than 2 second, you should talk to the team before committing (so we know to look into it for optimization)**

* The code should have been tested a bit, and please consider setting up a unit test.


## Installing the software

More information, installation instructions as well as usage tutorial, please see the wiki page. 
https://bitbucket.org/mpaterno/y3_cluster_cpp/wiki/browse/

* **There is a native macOS installation available, if you're willing to use a `devel` version of CosmoSIS. See [[Native macOS build]].**

* **Alternatively, you can install a docker version with instructions here https://bitbucket.org/mpaterno/y3_cluster_cpp/src/master/**

The docker installation is a pre-requisite -- unfortunately this is a commercially/community offered software, so we may not be able to help as much -- However, we expect the installation to be trivial for Mac and Ubuntu. If you use other operating systems, you may consider using a Ubuntu virtual machine for installation.
Note that docker more or less can only be installed on a machine that you have root access to (like your laptop), and the docker images can't be ran from an external hard disk.
After installing docker, there are several repos to clone and install, please follow the instructions on the bitbucket instruction page. This takes some downloading time but should be relatively straightforward.

* If you can't get the installation to work, email des-cosmosis-y3@fnal.gov.


-------------------------------------------
## Getting started

We devleop and test the code in the docker environment -- so the first thing would be to installthe docker software. 
The free version is sufficient for this purpose, and the installation instruction can be found here
https://www.docker.com/get-started

We will be installing four different `git` repositories, nested like Russian dolls. Stay on your toes!
_This_ repository is the very last one to clone.

### Get cosmosis-docker

In whatever directory is to become your top-level installation, clone the `cosmosis-docker` repository:

```bash
$ git clone https://<your-username>@bitbucket.org/mpaterno/cosmosis-docker.git
```

Make sure you are cloning `mpaterno/cosmosis-docker.git` and not `joezuntz/cosmosis-docker.git`. We have
some additions that are not (yet) in Joe's repository.

This will get some scripts, a `Dockerfile`, etc.

### Build your personal Docker image

Go into the git folder that you've just cloned.
The first build is:

```bash
$ ./get-cosmosis-and-vm cosmosis
```

If you get a failure when trying to download the `cosmosis-standard-library`, you will need to set up ssh key for your bitbucket account. Follow the instructions [here](https://confluence.atlassian.com/bitbucket/set-up-an-ssh-key-728138079.html#SetupanSSHkey-ssh2) 


This will clone the `cosmosis` and `cosmosis-standard-library` repositories, and will "run" your `Dockerfile` to build the docker image.

Please note that an important difference between `paterno/cosmosis-docker` and `joezuntz/cosmosis-docker` is
that `paterno/cosmosis-docker` checks our the correct version of the `cosmosis` and `cosmosis-standard-library` repositories.

### Get cuba_cpp

Move to the top-level `cosmosis` directory  and clone the repository:

```bash
$ git clone https://mpaterno@bitbucket.org/mpaterno/cubacpp.git
```

### Get y3_cluster_cpp

Move to the `cosmosis/cosmosis-standard-libary` directory, and in that directory clone our repository and likelihood module:

```bash
$ cd cosmosis/cosmosis-standard-library
$ git clone https://<your-username>@bitbucket.org/mpaterno/y3_cluster_cpp.git
$ git clone https://<your-username>@bitbucket.org/cosmosisclustercpp/cluster_mass_nc_like.git
```

Now you're ready to start building the software.

### Build cosmosis

From your top-level directory, start up the "VM" (really the Docker container):

```bash
$ ./start-cosmosis-vm cosmosis
```

This will land you at a _bash_ prompt, where you can build and then run CosmoSIS.
To build, just run _make_.

```bash
$ make
```

### Build the tests for our integration routines

Inside the running container, move to the `y3_cluster_cpp` directory, and
then build. Building is a two-step process: running CMake to generate many
makefiles, and then running _make_ to build the code.

```bash
$ cd cosmosis-standard-library/y3_cluster_cpp
$ cmake -DCMAKE_MODULE_PATH=/cosmosis/cubacpp/cmake/modules  -DCUBA_DIR=/usr/local -DCMAKE_BUILD_TYPE=Release .
$ make
$ ctest
```

The output from `ctest` should show some number of tests having run, and all passing.

## Keeping up-to-date

The `y3_cluster_cpp` project has many moving parts, and it is important to keep them synchronized.
As long as our code is still under rapid development, this will be challenging.
The following are suggestions on how to manage your use of the various repositories:

1. Keep your `cosmosis-docker` updated to the head of the master branch. This will help keep you using the correct versions of the Docker image layers.
2. Update `cosmosis` when you receive an email telling you to do so. We will update this infrequently, but it is important we are all using the same version.
3. Update `cosmosis-standard-library` regularly. We are using our own fork of Joe's repository so that we can both stay in sync with each other, and move forward more rapidly than does the official CosmoSIS release. We are working on the `develop` branch of this repository.
4. Update `y3_cluster_cpp` regularly. You should stay current with the head of the `master` branch. Unless you have a strong reason to do so, don't use other branches.

The current version of `cosmosis` is: 1.6rc1.

