package pi

import (
	"fmt"
	"math/rand"
	"runtime"
	"sync"
	"time"
)

func parallelExecution(
	loopUpperLimit int,
	worker func(r *rand.Rand, start, end int) float64,
) float64 {
	numCPU := runtime.NumCPU()
	base := loopUpperLimit / numCPU
	rem := loopUpperLimit % numCPU

	results := make([]float64, numCPU)
	var wg sync.WaitGroup

	for id := 0; id < numCPU; id++ {
		start := id * base
		end := start + base
		if id == numCPU-1 {
			end += rem
		}
		src := rand.NewSource(time.Now().UnixNano() + int64(id))
		r := rand.New(src)
		wg.Add(1)
		go workerLaunch(id, start, end, results, &wg, r, worker)
	}

	wg.Wait()
	total := 0.0
	for _, v := range results {
		total += v
	}
	return total
}

func workerLaunch(
	id, start, end int,
	results []float64,
	wg *sync.WaitGroup,
	r *rand.Rand,
	worker func(r *rand.Rand, start, end int) float64,
) {
	defer wg.Done()
	results[id] = worker(r, start, end)
}

func PerformPiCalculationRuns() {
	nRuns := 10
	runs := make([]float64, nRuns)

	for run := 0; run < nRuns; run++ {
		const nTests = 10_000_000_000
		startTime := time.Now()

		worker := func(r *rand.Rand, lo, hi int) float64 {
			count := 0
			for i := lo; i < hi; i++ {
				x := r.Float64()
				y := r.Float64()
				if x*x+y*y < 1.0 {
					count++
				}
			}
			return float64(count)
		}

		totalInCircle := parallelExecution(nTests, worker)
		approxPi := 4.0 * (totalInCircle / float64(nTests))
		elapsed := time.Since(startTime).Seconds()

		fmt.Printf("%8d %8.8f %8.8f\n", run, approxPi, elapsed)
		runs[run] = elapsed
	}

	sum := 0.0
	for _, v := range runs {
		sum += v
	}
	fmt.Printf("%8d %8.8f\n\n\n", nRuns, sum/float64(nRuns))
}
