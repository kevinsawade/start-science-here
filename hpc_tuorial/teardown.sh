#/bin/bash

docker-compose stop
docker-compose rm -f
docker volume rm hpc_tuorial_shared-vol
docker network rm hpc_tuorial_default
