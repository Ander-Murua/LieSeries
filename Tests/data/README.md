Put here the table files from https://www.ehu.eus/ccwmuura/bch.html:

    bchHall20.dat  bchLyndon20.dat  sbchHall19.dat  sbchLyndon19.dat

They are used by the test "Published-BCH-tables" in ../LieSeries.wlt, which reads
only the rows up to degree 8 (the first 71 rows of each file), so truncated copies
are enough. For example, in Mathematica:

    Do[Export[FileNameJoin[{"Tests", "data", f}],
         Take[Import["https://www.ehu.eus/ccwmuura/research/" <> f, "Table"], 71], "Table"],
       {f, {"bchHall20.dat", "bchLyndon20.dat", "sbchHall19.dat", "sbchLyndon19.dat"}}]
