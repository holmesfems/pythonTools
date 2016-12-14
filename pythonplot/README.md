python's plotting tools

usage:

    If you want to plot data for x-y1,y2...yn at one chart,just use plotfile_mal.py.

command:

    python plotfile_mal.py FILENAME OPTION_key=VALUE

example:

    python plotfile_mal.py plot1.list plot2.list

format of datafile:

    Line begin with '!' is the title line,and will be the title of each variable.

    Note that if the files have different name of variable X,the final name will be that
    occured at first file.

    Line begin with '#' is the comment line and will ignore this.

    Note that there is a option named commentInv that if you make this to true,the line begin
    with '#' will be recognized as data and others will be comment.(Inversion of the recognization
    of comment and data)

    Seperating columns can use '\t' or space

    The first column is data of X,and the other column is data of Y,for example:
    the row:
        1.0 2.0 3.0
    means [x,y1,y2] is [1.0,2.0,3.0]

    Value that mult to file means amplification of this data
