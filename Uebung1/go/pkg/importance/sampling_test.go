package importance

import (
	"math/rand"
	"testing"
)

func BenchmarkImportanceSampling(b *testing.B) {
	// choose a modest sample count so a single run is inexpensive—
	// we’ll let the harness call it b.N times.
	const nSamples = 100000

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

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		total := ParallelExecution(nSamples, worker)
		_ = total / float64(nSamples)
	}
}
