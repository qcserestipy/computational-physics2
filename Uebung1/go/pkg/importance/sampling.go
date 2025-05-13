package importance

import (
	"fmt"
	"math"
	"math/rand"
	"runtime"
	"sync"
)

func workerLaunch(id, lo, hi int, results []float64, wg *sync.WaitGroup, worker func(start, end int) float64) {
	defer wg.Done()
	results[id] = worker(lo, hi)
}

func ParallelExecution(
	nSamples int,
	worker func(start, end int) float64,
) float64 {
	numCPU := runtime.NumCPU()
	base := nSamples / numCPU
	rem := nSamples % numCPU

	var wg sync.WaitGroup
	results := make([]float64, numCPU)

	for id := 0; id < numCPU; id++ {
		start := id * base
		end := start + base
		if id == numCPU-1 {
			end += rem
		}
		wg.Add(1)
		go workerLaunch(id, start, end, results, &wg, worker)
	}
	wg.Wait()

	total := 0.0
	for _, partial := range results {
		total += partial
	}
	return total
}

func targetFunction(x float64) float64 {
	return math.Pow(x, -1.0/3.0) + x/10.0
}

func xFunction(u float64) float64 {
	return math.Pow((u+0.2746)/1.2746, 3.0/2.0)
}

func wFunction(x float64) float64 {
	return 0.8497 * math.Pow(x, -1.0/3.0)
}

func ImportanceSampling() {
	nSamples := 1000000

	worker := func(start, end int) float64 {
		sum := 0.0
		for i := start; i < end; i++ {
			u := rand.Float64()
			x := xFunction(u)
			w := wFunction(x)
			f := targetFunction(x)
			sum += f / w
		}
		return sum
	}

	grandTotal := ParallelExecution(nSamples, worker)
	estimate := grandTotal / float64(nSamples)
	fmt.Printf("nSamples = %d,  Ĩ = %f\n", nSamples, estimate)
}
