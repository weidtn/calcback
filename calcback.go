/*
Calculate the refractive index r_i = n + ik from the ellipsometric parameters Delta and Psi of a layer (300nm SiO2) on top of a substrate (Si-Wafer).
Input must be a .csv file
output will be a .png file (TODO and a .csv later?)

Currently the rerange and imrange are set for SiO2, which has 0 absorption (k), and n should be around 1.3-1.4
TODO: After I get the correct results, add a worker pool and goroutines to calculate every wavelenght in an own goroutine
*/
package main

import (
	"fmt"
	"io/ioutil"
	"log"
	"math"
	"math/cmplx"
	"strings"
	"time"

	"github.com/kniren/gota/dataframe"
	"github.com/kniren/gota/series"
	"gonum.org/v1/gonum/floats"
	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/plotutil"
	"gonum.org/v1/plot/vg"
)

const (
	phi_i           = 70 * math.Pi / 180 // incident angle in radians
	d_L             = complex(300,0)                // layer thickness [nm]
	n_air           = complex(1,0)                  // refractive index of air
	rerange         = 10               // real part from 0.1 to ...
	imrange         = 10                 // imaginary part from 0.1 to ...
	cmplx_pi        = complex(math.Pi, 0)
)

func timeTrack(start time.Time, name string) {
	elapsed := time.Since(start)
	log.Printf("%s took %s", name, elapsed)
}

func read_csv(file string) dataframe.DataFrame {
	defer timeTrack(time.Now(), "read_csv")
	content, err := ioutil.ReadFile(file)
	if err != nil {
		fmt.Print(err)
	}
	ioContent := strings.NewReader(string(content))
	df := dataframe.ReadCSV(ioContent,
		dataframe.WithDelimiter(','),
		dataframe.HasHeader(true))
	data := df.Select([]int{0, 1, 2, 5, 6}) // We only need the first three columns
	return data
}

func deltapsiplot(df dataframe.DataFrame) {
	psis := df.Col("Psi").Float()
	deltas := df.Col("Delta").Float()
	psi := make(plotter.XYs, df.Nrow())
	delta := make(plotter.XYs, df.Nrow())

	for i, floatlambda := range df.Col("lambda").Float() {
		psi[i].X = floatlambda
		psi[i].Y = psis[i]
		delta[i].X = floatlambda
		delta[i].Y = deltas[i]
	}
	p, err := plot.New()
	if err != nil {
		fmt.Println(err)
	}
	err = plotutil.AddLines(p, "Delta", delta, "Psi", psi)
	if err := p.Save(4*vg.Inch, 4*vg.Inch, "test.png"); err != nil {
		log.Fatal(err)
	}
	fmt.Println("Plot saved as test.png")
}

func Calc_rho(lambdafloat float64, n_S complex128) (n_rho [][]complex128) {
	lambda := complex(lambdafloat,0)
	var rho_L complex128
	var output [][]complex128
	// make a slice containing every possible n_L = n+ik
	// TODO: Put this outside of the function
	n := make([]float64, 100)
	k := make([]float64, 100)
	var nslice []complex128
	for _, x := range floats.Span(n, 1.0, rerange) {
		for _, y := range floats.Span(k, 0.001, imrange) {
			c := complex(x, y)
			nslice = append(nslice, c)
		}
	}
	//calculate for every n_L in nslice
	for _, n_L := range nslice {
		// Calculate Delta and Psi for given lambda

		// Snells law:
		phi_L := cmplx.Asin((cmplx.Sin(phi_i) * n_air) / n_L)
		phi_S := cmplx.Asin((cmplx.Sin(phi_L) * n_L) / n_S)

		// Fresnel equations:
		//
		// air/layer:
		rs_al := (n_air*cmplx.Cos(phi_i) - n_L  *cmplx.Cos(phi_L))/(n_air* cmplx.Cos(phi_i) + n_L  *cmplx.Cos(phi_L))
		rp_al := (n_L  *cmplx.Cos(phi_i) - n_air*cmplx.Cos(phi_L))/(n_L  * cmplx.Cos(phi_i) + n_air*cmplx.Cos(phi_L))

		// layer/substrate:
		rs_ls := (n_L*cmplx.Cos(phi_L) - n_S*cmplx.Cos(phi_S))/(n_L*cmplx.Cos(phi_L) + n_S*cmplx.Cos(phi_S))
		rp_ls := (n_S*cmplx.Cos(phi_L) - n_L*cmplx.Cos(phi_S))/(n_S*cmplx.Cos(phi_L) + n_L*cmplx.Cos(phi_S))

		beta := (2 * math.Pi / lambda) * d_L * n_L * cmplx.Cos(phi_L)

		rp_L := (rp_al + rp_ls * cmplx.Exp(1i*2*beta)) / (1 + rp_al*rp_ls*cmplx.Exp(1i*2*beta))

		rs_L := (rs_al + rs_ls * cmplx.Exp(1i*2*beta)) / (1 + rs_al*rs_ls*cmplx.Exp(1i*2*beta))

		rho_L = rp_L / rs_L
		row := []complex128{n_L, rho_L}
		output = append(output, row)
	}
	return output
}

func compare(n_rho [][]complex128, psi float64, delta float64) (n complex128) {
	// TODO: Check if this works as intended
	rho_giv := cmplx.Tan(complex(psi, 0)) * cmplx.Exp(complex(0,delta))
	var deltas []float64 // delta = difference between given and calculated rho
	for i := range n_rho {
		delta := cmplx.Abs(n_rho[i][1] - rho_giv)
		deltas = append(deltas, delta)
	}
	idx, _ := ArgMin(deltas)
	return n_rho[idx][0]
}

func ArgMin(array []float64) (int, float64) {
	// finds minimum in an array and returns index and value
	var min float64 = array[0]
	var idx int
	for i, value := range array {
		if value < min {
			min = value
			idx = i
		}
	}
	return idx, min

}

func plot_nk(lambdas, ns, ks []float64) {
	defer timeTrack(time.Now(), "plot_nk()")
	n := make(plotter.XYs, len(ns))
	k := make(plotter.XYs, len(ks))

	for i, lambda := range lambdas {
		n[i].X = lambda
		n[i].Y = ns[i]
		k[i].X = lambda
		k[i].Y = ks[i]
	}
	p, err := plot.New()
	if err != nil {
		fmt.Println(err)
	}
	err = plotutil.AddLines(p, "n", n, "k", k)
	if err := p.Save(4*vg.Inch, 4*vg.Inch, "nk.png"); err != nil {
		log.Fatal(err)
	}
	fmt.Println("Plot saved as nk.png")

}
func main() {
	defer timeTrack(time.Now(), "main")
	file := "300nmSiO2.csv"
	df := read_csv(file)
	nseries := series.New([]float64{}, series.Float, "n")
	kseries := series.New([]float64{}, series.Float, "k")
	//uncomment this and comment the for-loop to just calculate for 1 value
	// i := 0
	// lambda := 300.0
	for i, lambda := range df.Col("lambda").Float() {
		n_S := complex(df.Elem(i,3).Float(),df.Elem(i,4).Float()) // get refractive index of Substrate
		rhodata := Calc_rho(lambda, n_S)
		psi := df.Elem(i, 1).Float()
		delta := df.Elem(i, 2).Float()
		n := compare(rhodata, psi, delta)
		nseries.Append(real(n))
		kseries.Append(imag(n))
	}
	df = df.Mutate(nseries)
	df = df.Mutate(kseries)
	plot_nk(df.Col("lambda").Float(), df.Col("n").Float(), df.Col("k").Float())
	fmt.Println(df)
}
