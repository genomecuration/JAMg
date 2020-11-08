#Monitoring the Progress of Trinity

Since Trinity can easily take several days to complete a large trancriptome assembly job, it is useful to be able to monitor the process and to know at which stage (Inchworm, Chrysalis, Butterfly) Trinity is currently at.  There are a few general ways to do this:

- by running 'top', you'll be able to see which Trinity process is running and how much memory is being consumed.
- other downstream process will generate standard output.  Be sure to capture 'stdout' and 'stderr' when you run the Trinity script.  The format for capturing both stdout and stderr depends on your SHELL.  Figure out what shell you have by running:

      env | grep SHELL

    Using tcsh:

         Trinity ... opts ... > & run.log &

    Using bash:

        Trinity ... opts ... > run.log 2>&1 &

Note, under bash, to prevent the background process from being terminated once you close the shell, type 'exit' to leave the shell, or explore alternatives such as [nohup, disown, or screen](http://www.serverwatch.com/tutorials/article.php/3935306/Detach-Processes-With-Disown-and-Nohup.htm).

You can then 'tail -f run.log' to follow the progress of the Trinity throughout the various stages.
