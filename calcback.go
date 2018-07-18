package main

import (
	"fmt"
	"github.com/kniren/gota/dataframe"
	// "github.com/kniren/gota/series"
	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	// "gonum.org/v1/plot/plotutil"
	"gonum.org/v1/plot/vg"
	"io/ioutil"
	"log"
	"math"
	"math/cmplx"
	"strings"
	"time"
)

const phi_i = 70 * math.Pi / 180 // incident angle in radians
const d_L = 100                  // layer thickness [nm]
const step = 100                 // step size
const n_air float64 = 1          // refractive index of air
const n_S = 3.6449               // refractive index of substrate

func read_csv(file string) dataframe.DataFrame {
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

func save_plot(df dataframe.DataFrame) {
	// lambda := df.Col("Lambda").Int()
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
	p.Add(plotter.NewGrid())
	l1, _ := plotter.NewLine(psi)
	l2, _ := plotter.NewLine(delta)
	p.Add(l1, l2)
	if err := p.Save(4*vg.Inch, 4*vg.Inch, "test.png"); err != nil {
		log.Fatal(err)
	}
}

func calc_rho(data dataframe.DataFrame) (rho complex128, difference float64) {
	// Get variables from dataframe:
	lambda := data.Elem(0, 0).Float()
	delta := data.Elem(0, 1).Float()
	psi := data.Elem(0, 2).Float()
	var start float64 = 1.501
	// step := 0.001
	n_L := start
	rho_giv := cmplx.Tan(complex(psi, 0)) * cmplx.Exp(complex(0, delta))
	// total reflection:
	x := (math.Sin(phi_i) * n_air / n_L)
	if x > 1 || x < -1 || x == 0 {
		// continue
		return cmplx.NaN(), math.NaN()
	}
	// Calculate Delta and Psi for given lambda

	// Snells law:
	phi_L := math.Asin((math.Sin(phi_i) * n_air) / n_L)
	phi_S := math.Asin((math.Sin(phi_i) * n_air) / n_S)

	// Fresnel equations:

	// air/layer:
	rs_al := (n_air*math.Cos(phi_i) - n_L*math.Cos(phi_L)) / n_air * math.Cos(phi_i+n_L*math.Cos(phi_L))
	rp_al := (n_L*math.Cos(phi_i)-n_air*math.Cos(phi_L))/n_L*math.Cos(phi_i) + n_air*math.Cos(phi_L)

	// layer/substrate:
	rs_ls := n_L*math.Cos(phi_L) - n_S*math.Cos(phi_S)/n_L*math.Cos(phi_L) + n_S*math.Cos(phi_S)
	rp_ls := n_S*math.Cos(phi_L) - n_L*math.Cos(phi_S)/n_S*math.Cos(phi_L) + n_L*math.Cos(phi_S)

	beta := (2 * math.Pi / lambda) * d_L * n_L * math.Cos(phi_L)

	rp_L := (complex(rp_al, 0) * complex(rp_ls, 0) * cmplx.Exp(complex(0, 2*beta))) / (1 + complex(rp_al, 0)*complex(rp_ls, 0)*cmplx.Exp(complex(0, 2*beta)))

	rs_L := (complex(rs_al, 0) * complex(rs_ls, 0) * cmplx.Exp(complex(0, 2*beta))) / (1 + complex(rs_al, 0)*complex(rs_ls, 0)*cmplx.Exp(complex(0, 2*beta)))

	rho_L := rp_L / rs_L
	diff := cmplx.Abs(rho_L - rho_giv)
	fmt.Println(lambda)
	return rho_L, diff
}
func main() {
	file := "/home/aramus/Forschungspraktikum/tmm/300nmSiO2.csv"
	starttime := time.Now()
	df := read_csv(file)
	// fmt.Println(df)
	i := 0
	// for i := range df.Col("lambda").Float() {
	// fmt.Printf("%v  \n", calc_rho(df.Subset(i)))
	// }
	rho, diff := calc_rho(df.Subset(i))
	fmt.Printf("%.5f with a difference of %.5f \n", rho, diff)
	elapsed := time.Since(starttime)
	fmt.Printf("Done in %s", elapsed)

}
