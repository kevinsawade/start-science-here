 .. _install-docker-label:

 .. role:: bolditalic
  :class: bolditalic

.. role:: boldcode
  :class: boldcode

.. role:: italiccode
  :class: italiccode

====================================
Installation instructions for Docker
====================================

Set up the repository
---------------------

1. Update the ``apt`` package index and install packages to allow ``apt`` to use a repository over HTTPS:

.. code-block:: bash

   $ sudo apt-get update

   $ sudo apt-get install \
    apt-transport-https \
    ca-certificates \
    curl \
    gnupg \
    lsb-release

2. Add Docker’s official GPG key:

.. code-block:: bash

   $  curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /usr/share/keyrings/docker-archive-keyring.gpg

3. Use the following command to set up the **stable** repository. To add the **nightly** or **test** repository, add the word ``nightly`` or ``test`` (or both) after the word ``stable`` in the commands below. `Learn about nightly and test channels <https://docs.docker.com/engine/install/>`_.

.. note::

   The ``lsb_release -cs`` sub-command below returns the name of your Ubuntu distribution, such as ``xenial``. Sometimes, in a distribution like Linux Mint, you might need to change ``$(lsb_release -cs)`` to your parent Ubuntu distribution. For example, if you are using Linux Mint Tessa, you could use bionic. Docker does not offer any guarantees on untested and unsupported Ubuntu distributions.

.. code-block:: bash

   $ echo \
    "deb [arch=amd64 signed-by=/usr/share/keyrings/docker-archive-keyring.gpg] https://download.docker.com/linux/ubuntu \
    $(lsb_release -cs) stable" | sudo tee /etc/apt/sources.list.d/docker.list > /dev/null

Install Docker Engine
---------------------

1. Update the ``apt`` package index, and install the *latest version* of Docker Engine and containerd, or go to the next step to install a specific version:

.. code-block:: bash

   $ sudo apt-get update
   $ sudo apt-get install docker-ce docker-ce-cli containerd.io

2. Verify that Docker Engine is installed correctly by running the ``hello-world`` image.


.. code-block:: bash

   $ sudo docker run hello-world
