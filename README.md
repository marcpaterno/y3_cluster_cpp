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

