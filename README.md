This is the y3_cluster_cpp repository.

## Getting started

We will be installing four different `git` repositories, nested like Russian dolls. Stay on your toes!
_This_ repository is the very last one to clone.

### Get cosmosis-docker

In whatever directory is to become your top-level installation, clone the `cosmosis-docker` repository:

```bash
$ git clone https://<your-username>@bitbucket.org/mpaterno/cosmosis-docker.git
```

This will get some scripts, a `Dockerfile`, etc.

### Build your personal Docker image

The first build is:

```bash
$ ./get-cosmosis-and-vm cosmosis
```

This will clone the `cosmosis` and `cosmosis-standard-library` repositories, and will "run" your `Dockerfile` to build the docker image.

Now get the right branches of the `cosmosis` and `cosmosis-standard-library` repositories:

```bash
$ cd cosmosis
$ git checkout v1.5rc1
$ cd cosmosis-standard-library
$ get checkout v1.5rc1
```

### Get cuba_cpp

Move to the top-level `cosmosis` directory  and clone the repository:

```bash
$ git clone https://mpaterno@bitbucket.org/mpaterno/cubacpp.git
```

### Get y3_cluster_cpp

Move to the `cosmosis/cosmosis-standard-libary` directory, and in that directory clone our repository:

```bash
$ cd cosmosis/cosmosis-standard-library
$ git clone https://<your-username>@bitbucket.org/mpaterno/y3_cluster_cpp.git
```

Now you're ready to start building the software.
