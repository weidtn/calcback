package main

import (
	"fmt"
	"github.com/kniren/gota/dataframe"
	"github.com/kniren/gota/series"
	"gonum.org/v1/gonum/floats"
	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/plotutil"
	"gonum.org/v1/plot/vg"
	"io/ioutil"
	"log"
	"math"
	"math/cmplx"
	"strings"
	"time"
)

const (
	phi_i           = 70 * math.Pi / 180 // incident angle in radians
	d_L             = 100                // layer thickness [nm]
	n_air   float64 = 1                  // refractive index of air
	n_S             = 3.6449             // refractive index of substrate
	rerange         = 5                  // real part from 0.1 to ...
	imrange         = 25                 // imaginary part from 0.1 to ...
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
	data := df.Select([]int{0, 1, 2}) // We only need the first three columns
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

func calc_rho(lambda float64) (n_rho [][]complex128) {
	var rho_L complex128
	var output [][]complex128
	// make a slice containing every possible n_L = n+ik
	n := make([]float64, 100)
	k := make([]float64, 100)
	var nslice []complex128
	for _, x := range floats.Span(n, 0.1, rerange) {
		for _, y := range floats.Span(k, 0.1, imrange) {
			c := complex(x, y)
			nslice = append(nslice, c)
		}
	}

	//calculate for every n_L in nslice
	for _, n := range nslice {
		n_L := real(n)
		// n_L := n
		// TODO total reflection:
		// Use real only for Snells law?
		//
		x := (math.Sin(phi_i) * n_air / n_L)
		if x > 1 || x < -1 || x == 0 {
			continue
			//  append NaN?
		}
		// Calculate Delta and Psi for given lambda

		// Snells law:
		phi_L := math.Asin((math.Sin(phi_i) * n_air) / n_L)
		phi_S := math.Asin((math.Sin(phi_i) * n_air) / n_S)

		// Fresnel equations:
		//
		// air/layer:
		rs_al := (n_air*math.Cos(phi_i) - n_L*math.Cos(phi_L)) / n_air * math.Cos(phi_i+n_L*math.Cos(phi_L))
		rp_al := (n_L*math.Cos(phi_i)-n_air*math.Cos(phi_L))/n_L*math.Cos(phi_i) + n_air*math.Cos(phi_L)

		// layer/substrate:
		rs_ls := n_L*math.Cos(phi_L) - n_S*math.Cos(phi_S)/n_L*math.Cos(phi_L) + n_S*math.Cos(phi_S)
		rp_ls := n_S*math.Cos(phi_L) - n_L*math.Cos(phi_S)/n_S*math.Cos(phi_L) + n_L*math.Cos(phi_S)

		beta := (2 * math.Pi / lambda) * d_L * n_L * math.Cos(phi_L)

		rp_L := (complex(rp_al, 0) * complex(rp_ls, 0) * cmplx.Exp(complex(0, 2*beta))) / (1 + complex(rp_al, 0)*complex(rp_ls, 0)*cmplx.Exp(complex(0, 2*beta)))

		rs_L := (complex(rs_al, 0) * complex(rs_ls, 0) * cmplx.Exp(complex(0, 2*beta))) / (1 + complex(rs_al, 0)*complex(rs_ls, 0)*cmplx.Exp(complex(0, 2*beta)))

		rho_L = rp_L / rs_L
		row := []complex128{n, rho_L}
		output = append(output, row)
	}
	return output
}

func compare(n_rho [][]complex128, psi float64, delta float64) (n complex128) { // TODO Complex abs is not real abs.
	// n_rho [][0] contains n_L
	// n_rho [][1] contains rho
	rho_giv := cmplx.Tan(complex(psi, 0)) * cmplx.Exp(complex(0, delta))
	var deltas []float64 // delta = difference between given and calculated rho
	for i := range n_rho {
		delta := cmplx.Abs(n_rho[i][1] - rho_giv) //Problem is here
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
	// fmt.Println(df)
	// deltapsiplot(df)
	nseries := series.New([]float64{}, series.Float, "n")
	kseries := series.New([]float64{}, series.Float, "k")
	// i := 0
	// lambda := 300.0
	for i, lambda := range df.Col("lambda").Float() {
		rhodata := calc_rho(lambda)
		delta := df.Elem(i, 1).Float()
		psi := df.Elem(i, 2).Float()
		n := compare(rhodata, psi, delta)
		// fmt.Println(rhodata)
		nseries.Append(real(n))
		kseries.Append(imag(n))
		// fmt.Println(df.Elem(i, 0), real(n))
		// fmt.Println(real(n))
		// fmt.Println(imag(n))
	}
	df = df.Mutate(nseries)
	df = df.Mutate(kseries)
	plot_nk(df.Col("lambda").Float(), df.Col("n").Float(), df.Col("k").Float())
	fmt.Println(df)
}
