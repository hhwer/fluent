
set term post eps color solid 
set output "fig1.eps"
set size 1,1
set multiplot layout 2,2
set size 0.5,0.5

a=4

U = "norm001U".a.".txt"
W = "norm001W".a.".txt"
set origin 0,0
set title "tau=0.01,N=".(2*a)
set xlabel 't(s)'
set ylabel 'norm_inf'
set xtics("0"0,"100"10000,"200"20000,"300"30000,"400"40000,"500"50000,"600"60000,"700"70000,"800"80000,"900"90000,"1000"100000)
plot U with line lt 1  , W with line lt -1

set origin 0,0.5 
set size 0.5,0.5
set xtics("0"0,"100"1000,"200"2000,"300"3000,"400"4000,"500"5000,"600"6000,"700"7000,"800"8000,"900"9000,"1000"10000)
U = "norm01U".a.".txt"
W = "norm01W".a.".txt"
set title "tau=0.1,N=".(2*a)
set xlabel 't(s)'
set ylabel 'norm_inf'
plot U with line lt 0 lw 2 , W with line lt -1

a=4

U = "norm001U".a.".txt"
W = "norm001W".a.".txt"
set origin 0.5,0
set title "tau=0.01,N=".(2*a)
set xlabel 't(s)'
set ylabel 'norm_inf'
set xtics("0"0,"100"10000,"200"20000,"300"30000,"400"40000,"500"50000,"600"60000,"700"70000,"800"80000,"900"90000,"1000"100000)
plot U with line lt 0 lw 2 , W with line lt -1

set origin 0.5,0.5 
set xtics("0"0,"100"1000,"200"2000,"300"3000,"400"4000,"500"5000,"600"6000,"700"7000,"800"8000,"900"9000,"1000"10000)
U = "norm01U".a.".txt"
W = "norm01W".a.".txt"
set title "tau=0.1,N=".(2*a)
set xlabel 't(s)'
set ylabel 'norm_inf'
plot U with line lt 0 lw 2 , W with line lt -1

unset multiplot
set output

