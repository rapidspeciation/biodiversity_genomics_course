# Working on the cluster

## Login
The login process is in multiple steps.

### Step 1 – log in to the cluster

First we need to log in to the cluster with authentification:

```
kn9@mib118791s ~ % ssh kn9sanger_ac@168.176.34.122
```

The first time you log in you will be asked if you want to continue: Are you sure you want to continue connecting (yes/no/[fingerprint])?
Type yes to add the server to known hosts.

Then it will ask for your password, which will be hidden when you type or paste it in, so just add you password and press enter, then it will connect. 

```
The authenticity of host '168.176.34.122 (168.176.34.122)' can't be established.
ED25519 key fingerprint is: SHA256:Rki19uJ+LyynhL3VHP1eOG+mUHV5FMVqoFDhj9EcYUI
This key is not known by any other names.
Are you sure you want to continue connecting (yes/no/[fingerprint])? yes
Warning: Permanently added '168.176.34.122' (ED25519) to the list of known hosts.
** WARNING: connection is not using a post-quantum key exchange algorithm.
** This session may be vulnerable to "store now, decrypt later" attacks.
** The server may need to be upgraded. See https://openssh.com/pq.html
kn9sanger_ac@168.176.34.122's password: 
Linux cecc.unal.edu.co 5.10.0-45-amd64 #1 SMP Debian 5.10.259-1 (2026-07-02) x86_64
 
The programs included with the Debian GNU/Linux system are free software;
the exact distribution terms for each program are described in the
individual files in /usr/share/doc/*/copyright.
 
Debian GNU/Linux comes with ABSOLUTELY NO WARRANTY, to the extent
permitted by applicable law.
Last login: Wed Jul  8 10:48:32 2026 from 186.28.174.88
kn9sanger_ac@cecc:~$ 

```

The first step has taken us to the login server cecc

```
kn9sanger_ac@cecc:~$ ls
kn9sanger_ac@cecc:~$ pwd
/homes/C_Computacion/kn9sanger_ac
```

### Step 2

To do work and save files we need to use a server with more resources, for this we need a different server, called perseus.
From the login node we will connect to the perseus server

```
kn9sanger_ac@cecc:~$ ssh perseus
The authenticity of host 'perseus (168.176.8.19)' can't be established.
ECDSA key fingerprint is SHA256:1U6aFVnYu9wWeVAusNGaEu3DM0bbg8WVpXGrH87cj+o.
Are you sure you want to continue connecting (yes/no/[fingerprint])? yes
Warning: Permanently added 'perseus,168.176.8.19' (ECDSA) to the list of known hosts.
Linux perseus 5.10.0-30-amd64 #1 SMP Debian 5.10.218-1 (2024-06-01) x86_64
 
The programs included with the Debian GNU/Linux system are free software;
the exact distribution terms for each program are described in the
individual files in /usr/share/doc/*/copyright.
 
Debian GNU/Linux comes with ABSOLUTELY NO WARRANTY, to the extent
permitted by applicable law.
Last login: Wed Jul  8 10:48:36 2026 from 168.176.34.122
kn9sanger_ac@perseus:~$
```

Check where you are with pwd:

```
kn9sanger_ac@perseus:~$ pwd
/homes/C_Computacion/kn9sanger_ac

```

### Step 3 – computational resources

We do not have much computational power or space here in the login node of perseus.


Alternative 1:

Use hercules2 - 6 

This works for simple processes and we do not need to allocate space, we can just connect via ssh. 

```
kn9sanger_ac@perseus:~$ ssh hercules2

```


Alternative 2:

hercules are slow so for more computationally intensive tasks we should use a different node system

To list the nodes with more computational power 

```
kn9sanger_ac@perseus:~$
kn9sanger_ac@perseus:~$ sinfo -Mfisica

```

For these we first have to allocate space on a node:

```
kn9sanger_ac@perseus:/scratchsan/C_computacion/kn9sanger_ac$ salloc -Mfisica -pcpu.cecc
salloc: Pending job allocation 201342585
salloc: job 201342585 queued and waiting for resources
salloc: job 201342585 has been allocated resources
salloc: Granted job allocation 201342585
salloc: Waiting for resource configuration
salloc: Nodes boltzmann are ready for job


```

Default is cpu 2, 1 GB, but we can specify more by ?

To see what is allocated use squeue

```

kn9sanger_ac@perseus:/scratchsan/C_computacion/kn9sanger_ac$ squeue -u kn9sanger_ac
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
         201342585  cpu.cecc interact kn9sange  R       0:21      1 boltzmann


```


And if we need more resources, we have this group of nodes that we also can use
sinfo -Mcecc

To use the resources in the node we need to ssh to that machine, in this case boltzmann but it could be other names depending on the allocation.

```
kn9sanger_ac@perseus:/scratchsan/C_computacion/kn9sanger_ac$ ssh boltzmann

kn9sanger_ac@boltzmann:~$ pwd
/homes/C_Computacion/kn9sanger_ac

```

### Step 4 – move to working directory

We are using the node Boltzmann, which has the computational resources we need, but are located in our home directory, which is very small. We need to move to a directory with more space.

All the work should be done under the directory called scratchsan, it has 100GB per user. (not available from the login node only from perseus).

Move to the scratchan directory:

```
kn9sanger_ac@boltzmann:~$ cd /scratchsan/C_computacion/kn9sanger_ac/
kn9sanger_ac@boltzmann:/scratchsan/C_computacion/kn9sanger_ac$ pwd
/scratchsan/C_computacion/kn9sanger_ac

```

Here is where we can save files and folders and work with constructing scripts and run jobs.

Now we can finally start working.


## Software

To access software:

```
kn9sanger_ac@boltzmann:/scratchsan/C_computacion/kn9sanger_ac$ module load envs/anaconda3

(base) kn9sanger_ac@boltzmann:/scratchsan/C_computacion/kn9sanger_ac$ conda env list


conda activate <software>

conda deactivate

```

Check software and version


## Upload and download data

To upload a file or folder **from your local computer** to the cluster use the -J flag to jump across the login node to the university cluster:

scp -r -J user@servername file_to_upload user@servername2:absolute_path/to/directory/

You might have to add your password twice.

Make sure you are in your local computer!

```bash

scp -r -J kn9sanger_ac@168.176.34.122 test_up.txt kn9sanger_ac@perseus:/scratchsan/C_computacion/kn9sanger_ac/intro_unix/

```

To download to your computer:

scp -r -J user@servername user@servername2:absolute_path/to/directory/file_to_download /path/to/target/

Make sure you are in your local computer!

```bash

scp -r -J kn9sanger_ac@168.176.34.122 kn9sanger_ac@perseus:/scratchsan/C_computacion/kn9sanger_ac/intro_unix/testfile.txt ./

```





