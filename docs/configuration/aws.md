# nfcore/rnaseq: AWS Configuration

## Introduction

Amazon Web Services (AWS) are cloud-based compute systems, suitable for large-scale bioinformatics analysis. They can be an attractive way to run scalable analyses without needing to pay large up-front expenses for hardware and maintenance.

nfcore/rnaseq is built with NextFlow, which has a number of AWS integrations built in. However, AWS is vast and almost infinitely configurable, and as such there are several different ways in which you can use it for your analyses. Which approach you use depends on your needs - the frequency with which you run the pipeline, the scale of data you need to analysis and your familiarity with running in the cloud.

In this document, we describe the steps taken when running nfcore/rnaseq on AWS hardware using three different methods. These have been written as a walk-through based on our experience at the time of writing. Note that we have used AWS based in Ireland (`eu-west-1`) and that the availability of different hardware can vary according to location.

Please note that the state of Nextflow AWS integration has progressed very quickly, and some of these instructions may be out of date.

If you find any problems with this documentation, please let us know or better still, submit a pull-request with a fix.

> **Update - December 2017**
> Nextflow has recently added support for the new AWS Batch system. This method is probably better than those described above, but we haven't tried it yet so don't have any documentation written. For more information, see https://www.nextflow.io/blog/2017/scaling-with-aws-batch.html

## Table of Contents

1. [AWS Concepts](#aws-concepts)
    * Basic introduction to some AWS terms that you'll need to know.
2. [Fully Manual Single Node Setup](#1-fully-manual-single-node-setup)
    * The most basic way to run the pipeline. Involves creating a server on AWS and running NextFlow without any clever parallelisation stuff.
    * Good when you're just getting started, or if you have a small amount of data.
3. [NGI Automated Single Node Setup Script](#2-ngi-automated-single-node-setup-script) _(coming soon)_
4. [NextFlow Integration - Elastic Clusters](#3-nextflow-integration---elastic-clusters)

## AWS Concepts
The people behind AWS love acronyms. There is a huge amount of jargon to get your head around when starting with this. We've tried to keep this documentation as simple as possible, but it helps to understand a few basics:

### Instance Types
Broadly speaking, everything in this document will use compute resources that are part of the AWS EC2 - that is, the Amazon Web Services Elastic Cloud 2. These are on-demand servers that you hire on a use-by-time basis. Each one you run is called an _instance_. You can choose different instance types according to your requirements as each has different hardware specifications. See [Amazon EC2 Instance Types](https://aws.amazon.com/ec2/instance-types/) for a full list.

### Storage types
As well as an instance to run computation work on, we'll also need some file storage. This is organised separately from the instance and is then mounted to it. There are three types which are worth knowing about:

1. EFS _(Elastic File System)_
    * This is the most expensive filesystem. It can scales to whatever size of data you need, without any need to predefine a requirement. You just pay for what you use. A blank EFS filesystem is free to keep alive.
    * EFS filesystems can be mounted on multiple different instances.
2. EBS _(Elastic Block Store)_
    * This filesystem is cheaper than EFS. You can choose what type of storage to use (eg. solid-state or hard disk) for varying prices.
    * Unlike EFS, you need to know how much storage will be required in advance. If you go over this then everything will break.
    * EBS volumes can only be mounted on a single instance.
3. S3 _(Simple Storage Service)_
    * Cheap long term storage.
    * It's good for archiving data to once you're finished, but no good for running any analysis on.

In this documentation, we use EFS for running all of our analysis. There are a few benefits it gives - if the instance dies, we don't lose our data. It's faster and we don't need to worry about how much space we're going to need. The downside is that it's expensive, so you should clear it out as soon as you're done.

### Regions
AWS resources are based in data centres, which are based in regions. Typically, all data and hardware is restricted to a single region. These examples are written based on experience with Ireland. They should work in US-West with some tweaking. Other regions may have varying hardware available, so approach with care.

## 1. Fully Manual Single Node Setup
Simplest to setup, just using a single instance.

1. Log into your AWS console
    * Select the Ireland region (`eu-west-1`)
    * This can all be done either with the web interface or the command line client. Here we will describe the web interface.
2. Launch an EFS File system
    1. Go to the _Storage: EFS_ Dashboard
    2. Click _Create file system_
    3. Leave everything as default and click through the pages to create the file system.
        * Note the name of the security group that is created
    4. Once created, expand the EFS details and copy the _DNS name_ to your clipboard for later.
3. Set up a new security group
    1. Go to the EC2 Dashboard
    2. Select _Security Groups_ under the _Network & Security_ section in the left navigation
    3. Click _Create Security Group_
    4. Set the name and description to `SSH` (or whatever you like)
    5. Click _Add Rule_
    6. Select `SSH` as _Type_ and `Anywhere` as _Source_ (or a more limited IP range if you prefer)
    7. Click _Create_
4. Launch an Amazon Machine Instance (AMI)
    1. Go to the EC2 Dashboard
    2. Click _Launch Instance_
    3. Select the _Community AMIs_ tab and search for the nfcore/rnaseq AMI that we have created
        * AMI id: `ami-9a5c79e9`
        * This image contains all of the required software for the nfcore/rnaseq pipeline and is maintained by us for easy use by the community.
    4. Click Select
    5. Select the `m4.4xlarge` Instance Type
        * If running with STAR and Human data, we need around 34GB RAM. This is the cheapest instance type with that memory available.
        * If aligning with HISAT2 or using a smaller genome then you can use smaller instances, such as `m4.2xlarge` or `t2.2xlarge` (`t2` is less efficient but cheaper).
    6. Click _Configure Instance_
    7. Optional: Click _Request Spot Instances_
        * Enter a maximum bid value. Choosing something close to or above the on-demand price will mean that you rarely if ever get kicked off, but still get reduced costs most of the time.
        * Using spot pricing typically means that you save around 87% of the list price. See the [AWS Spot Bid Advisor](https://aws.amazon.com/ec2/spot/bid-advisor/) for more information.
        * For `m4.4xlarge`, $0.95 should be sufficient at time of writing.
    8. Select one of the subnets and note the name.
    9. Expand _Advanced Details_ at the bottom.
        * Under _User data_ add the following as _as text_:
        ```yaml
        # cloud-config
        package_upgrade: true
        packages:
        - nfs-utils
        runcmd:
        - mkdir -p /mnt/efs
        - chown ec2-user:ec2-user /mnt/efs
        - echo "SUBNET-ID.YOUR-DNS-NAME-HERE:/ /mnt/efs nfs4 nfsvers=4.1,rsize=1048576,wsize=1048576,hard,timeo=600,retrans=2 0 0" >> /etc/fstab
        - mount -a -t nfs4
        - sudo chown ec2-user /mnt/efs
        ```
        * Replace `SUBNET-ID` and `YOUR-DNS-NAME-HERE` with the subnet ID from above and the DNS name that you copied earlier.
        * For example, you may copy in something like: `eu-west-1c.fs-hwjer3hs.efs.eu-west-1.amazonaws.com`
    10. Click through the tabs until you get to _Configure Security Group_
    11. Use _Select an existing security group_ and choose the two that you created earlier
        * One to access the EFS, one to SSH into the instance.
    12. Click _Review and Launch_, then _Launch_
    13. Create a new Key Pair if you don't already have one
        * This is used to access the instance by `ssh`
    14. Click _Request Spot Instances_ :tada:
5. Access the new instance
    1. Once the spot instance has been given, go to _Instances_ on the EC2 dashboard.
    2. Copy the _IPv4 Public IP_ from the new instance details tab
    3. Load up a console on your local computer (eg. _Terminal_)
    4. Once the _Instance State_ and _Status Checks_ are complete, then SSH into the new instance
        * `ssh -i /path/to/key.pem ec2-user@YOUR-IP`
        * Replace `YOUR-IP` with the _IPv4 Public IP_ that you copied above.
    5. If you manage to log in successfully, pat yourself on the back and have a cookie. :cookie:
6. Get Nextflow running with the pipeline
    1. Go to the `efs` directory and work there: `cd /mnt/efs`
    2. Download the required references
    3. Download the input data
        * You can try this with our test data:
        ```
        curl https://export.uppmax.uu.se/b2013064/test-data/ngi-rna_test_set.tar.bz2 > ngi-rna_test_set.tar.bz2
        tar xvjf ngi-rna_test_set.tar.bz2
        ```
    4. Run the pipeline!
        * If you've downloaded the test data, this command should work:
        ```
        nextflow run nf-core/rnaseq -profile base --gtf ngi-rna_test_set/genes.gtf --bed12 ngi-rna_test_set/genes.bed --star_index ngi-rna_test_set/star/ --reads "ngi-rna_test_set/*.fastq.gz"
        ```
6. Retrieve your finished results to your local computer
    1. Use SCP on your local machine to pull your data back from the EFS file system:
    ```
    scp -r -i /path/to/key.pem ec2-user@YOUR-IP:/mnt/efs/results .
    ```
7. Stop everything and clean up
    1. Wipe the EFS filesystem so that it doesn't cost any money
        * On the AWS instance, delete all files on EFS:
        ```
        cd ~
        rm -rf /mnt/efs/*
        ```
        * Now it's empty, it should be free to keep it around (so no need to delete it)
    2. Go to the EC2 Dashboard
    3. Select the running instance and click _Actions_ > _Instance State_ > _Terminate_
        * **If you leave this running, things will get expensive really quickly**
    4. Cancel the spot request under _Spot Requests_



## 2. NGI Automated Single Node Setup Script
Denis wrote a magic script which essentially automates nearly all of the above to make
your life easier. More extensive docs and stuff coming soon.

## 3. NextFlow Integration - Elastic Clusters
Nextflow has built-in support for working with the AWS cloud. It can do many of the things
described above and has native support for things like S3, auto-scaling elastic clusters
and EFS mounts.

We have created a pipeline profile called `aws` (specify with `-profile aws`). This contains
a configuration designed for running on the AWS cloud, including links to common reference
genomes that we store for you on S3. This means that running on the cloud is very similar
to running on a HPC system and is super quick and easy!

### AWS Setup (one time only)
Unfortunately, we can't specify everything for you in the pipeline, you have to create
a few things on amazon first. These things only need to be done once and can be
reused for multiple pipeline runs.

#### Set up authentication
There are a few ways to do this, but we recommend setting up a config as used by the
AWS command line tools, which Nextflow also understands. See the
[AWS docs](http://docs.aws.amazon.com/cli/latest/userguide/cli-chap-getting-started.html#cli-config-files)
on this topic. For alternative methods to do this, see the
[Nextflow docs](https://www.nextflow.io/docs/latest/amazons3.html#aws-access-and-secret-keys).

We recommend that you use IAM on Amazon to limit the permissions for the user that
runs the pipelines. Ensure that this user has access to EC2 and S3 as a minimum.

#### Create an S3 bucket for your results
We recommend using S3 to store your finished pipeline results, and also the intermediate
files (nextflow's `work` directory). We can do a clever trick to tell S3 to automatically
delete the work directory files after a certain amount of time (eg. a week). This way,
you still have the intermediate files if you want to re-run the pipeline for any reason,
but they won't cost you lots of money if you forget to delete them when you're finished.

Go to your [S3 AWS console page](https://console.aws.amazon.com/s3/home) and create a new
bucket (if required). Go into the bucket and create a directory called `results` and a
directory called `work`. Back in the top level view, select the bucket and click the
_Properties_ tab. Expand the dropdown called _Lifecycle_ and then _Add Rule_. Set the
target as `work/` (or whatever you called your folder) and configure to
_Permanently Delete_ objects after a certain number of days (eg. `7`). It doesn't hurt
to also select _End and Clean up Incomplete Multipart Uploads_ whilst you're here.
Give the rule a name and hit save.

Finally, make a note of the bucket name and directories created. You'll need these later.

#### Create a VPC + Subnet
Nextflow requires a VPC (_virtual private cloud_) to run properly. You may have one that
has already been automatically created, or you may need to make a new one. Either way,
you need to make a note of the subnet ID.

Go to your [VPC dashboard page](https://eu-west-1.console.aws.amazon.com/vpc/home?region=eu-west-1)
(link is for `eu-west-1` - check that you're in the correct region). The dashboard will say
whether you have any subnets. If you do, click _Subnets_ in the navigation on the left
and make a note of one of the _Subnet ID_ values (_eg._ `subnet-05222a43`).

#### Ensure that your security group allows access
In order for you to log into your cluster, you need to tell AWS to allow you access.
To do this, go to the [Security Groups tab](https://eu-west-1.console.aws.amazon.com/ec2/v2/home?region=eu-west-1#SecurityGroups:sort=groupId)
under the EC2 dashboard and select the default security group. Switch to the
_Inbound_ tab and check that you have permissions to access resources. If in doubt,
set `All traffic` to `Anywhere` - not very secure, but should work (your cluster
will still be secured with your personal SSH key).

#### Configure Nextflow
We now need to configure Nextflow so that it knows how to create your cluster.
As it is personal to you, it cannot go in the main pipeline config. Instead, we recommend
saving everything in your personal config file, typically found in your home directory at
`~/.nextflow/config`. Open this file and add the following:

```groovy
cloud {
  imageId = 'ami-43f49030'
  instanceType = 'm4.large'
  subnetId = 'YOUR-SUBNET-ID'
  spotPrice = 1
  autoscale {
    enabled = true
    minInstances = 1
    maxInstances = 10
    instanceType = 'm4.2xlarge'
    spotPrice = 1
    terminateWhenIdle = true
  }
}
process {
    // Set here according to resources available on `m4.2xlarge`
    // Memory should be a little less than Amazon lists for node type
    cpus = 8
    memory = 31.GB
}
```

Remember to replace `YOUR-SUBNET-ID` with the subnet that you created / copied the address
for above.

#### A note about some pipeline defaults
Note that the above details are what we consider to be "sensible defaults". These settings
relate to the resources that you will use and how much that you'll be willing to pay for them.
Notable are `cloud.instanceType` and `cloud.autoscale.instanceType`, which define what type
of EC2 instance we're going to use for the head node and worker nodes. We've set `m4.large`
and `m4.2xlarge` respectively, which should work well for most Human runs. If your samples
are from an organism with a small genome, you may be able to get away with smaller and
cheaper instances.

Secondly, we define `cloud.spotPrice` and `cloud.autoscale.spotPrice`. These define our
"maximum bids" for the instances on the [AWS spot market](https://aws.amazon.com/ec2/spot/).
Both of these values are set to a super-high $1 per instance per hour. This means that the
chance of you being outbid by someone else is tiny, and that your instances are very unlikely
to be stopped. Note that it _doesn't_ mean that you're going to end up paying that much -
you will be charged the lowest possible spot market rate at the time that you run.
You can run the `nextflow cloud spot-prices` command to see the current spot price values.

### Upload your raw data to S3
In order to run your pipeline, Nextflow needs to be able to access your raw sequencing
data. The best way to do this is to create a directory in the same S3 bucket you created
above and to use the [AWS s3 command line tools](http://docs.aws.amazon.com/cli/latest/userguide/using-s3-commands.html)
to upload your data.

For example, if you created a directory inside the s3 bucket above called
`raw_data/my-project`, you could sync a local directory with the following command:

```bash
aws s3 sync my_data/ s3://my-bucket/raw_data/my_project/
```

### Create your cluster with nextflow
Now that everything is in place, we can get on with creating an AWS instance to use
for our pipeline run! Nextflow handles most of this for us (see these two blog posts
for more details: [one](https://www.nextflow.io/blog/2016/deploy-in-the-cloud-at-snap-of-a-finger.html)
and [two](https://www.nextflow.io/blog/2016/enabling-elastic-computing-nextflow.html)).

The above config should define basically all of the settings, so we just need to run
the following command:

```bash
nextflow cloud create MY-CLUSTER-NAME
```

You should see a confirmation with all of the config details that we saved above.
Type `y` and enter, then wait for your cluster to be started.

Once complete, copy the `ssh` command printed by Nextflow and run it to log into
your new cluster.

### Start your analysis run
Once you're logged into your cluster, you run the pipeline much as you normally
would, with a couple of additional parameters. For example:

```
~/nextflow run nf-core/rnaseq -w s3://my-bucket/work --reads 's3://my-bucket/raw_data/my-project/sample_*_{1,2}.fastq' --outdir 's3://my-bucket/results/my-project' -profile aws --genome GRCh37
```

A description of these parameters:

* `~/nextflow`
    * Your cluster comes with an installation of Nextflow. Make sure you use this path.
* `-w`
    * Sets the work directory for Nextflow. We tell it to use our s3 bucket work directory, where
    files are automatically deleted after one week.
    * NB: _one_ hyphen
* `--reads`
    * Our input data, that we've uploaded to s3 with the aws command line tools.
    * Remember to have in quotes, as always!
* `--outdir`
    * Our output directory, also on s3 but in our results directory where they're safe for long term storage.
* `-profile aws`
    * Tells the pipeline to use the `aws` config file.
    * NB: _one_ hyphen
* `--genome`
    * Key for one of the iGenome references. These are also kept on s3, but in a bucket that we host for you. These paths are specified in the `aws` profile config file.

### Shut everything down and get your results
When your pipeline has finished running, you need to shut your cluster down. Log out from the cloud
and then run the following command:

```bash
nextflow cloud shutdown MY-CLUSTER-NAME
```

All of your results are stored on s3, you can now sync these to your local computer if you'd like:

```bash
aws s3 sync s3://my-bucket/results/my_project/ results/
```

If you're happy with everything, you can also clean our the s3 work directory so that it doesn't
cost money by hosting the intermediate files for a week.

### Troubleshooting

#### Warnings about resources not being available
When you run your pipeline, you may get log messages that look like the following:

```
WARN: ### Task (id=6) requests an amount of resources not available in any node in the current cluster topology -- CPUs: 1; memory: 8 GB
WARN: ### Task (id=7) requests an amount of resources not available in any node in the current cluster topology -- CPUs: 2; memory: 16 GB
WARN: ### Task (id=8) exceed the number of CPUs provided by autoscaling instance type: m4.2xlarge -- req: 10; provided: 8
WARN: ### Task (id=9) exceed the number of CPUs provided by autoscaling instance type: m4.2xlarge -- req: 10; provided: 8
```

These are telling you that some of the processes in the pipeline are asking for impossible
numbers of cpus or memory. To fix this, you need to set the `params.cpus` and `params.memory`
in the config file described above to numbers that are available on your worker node type.
Note that the memory should be a little _below_ what Amazon lists as the available capacity.
