This is the y3_cluster_cpp repository.

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

This will clone the `cosmosis` and `cosmosis-standard-library` repositories, and will "run" your `Dockerfile` to build the docker image.

Now get the right branches of the `cosmosis` and `cosmosis-standard-library` repositories:

```bash
$ cd cosmosis
$ git checkout develop
$ cd cosmosis-standard-library
$ git checkout v1.5rc1
```

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
