@page cicd_dev CI/CD

@section new_deb Creating a new Debian repository

The repository must first be established with a tool such as 
[aptly](https://www.aptly.info/). The repo is hosted on `192.168.1.122` (aka `nemosys-repository.illinois.rocstar`). The user to create the repo is `administrator`, while the user that holds the files is a non-admin user `nemosys`. The commands to create, publish, and add to the repo would be similar to those shown below:

    aptly repo create -architectures="x86_64" -component="main" -distribution="focal" nemosys-repository-focal
    aptly publish repo -gpg-key=$KEY_ID -architectures="x86_64,i386,amd64" nemosys-repository-focal
    aptly repo add nemosys-repository-focal /home/nemosys/repo-focal

The last command updates the repository with debs in that directory.

A webserver reverse proxy must be set up to serve the location that has the repository, if it has not already been handled as in the bionic version. An example that used [nginx](https://www.nginx.com/) is present on the server.

The `administrator` user has a cron job that updates the repos when new files are added. The script is `/home/administrator/repo-updater.sh`, and it must be used to update the new repo.

The final step is to update the <em>NEMoSys</em> project `gitlab-ci.yml` file to deploy the correct debs to the `/home/nemosys/repo-focal` folder.