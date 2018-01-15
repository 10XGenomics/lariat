// Copyright (c) 2015 10X Genomics, Inc. All rights reserved.

package optimizer

import (
	"math"
	"math/rand"
)


type Optimizable interface {
	GenerateMove(accept_move func(log_p_curr float64, log_p_next float64) bool) Optimizable
}

func Optimize(current_model Optimizable, start_temp float64, temperature_steps int, steps_per_temp int) Optimizable {
	var seed int64 = 1
	var random = rand.New(rand.NewSource(seed))
	temp_steps := GetExponentialTemperatureSteps(start_temp, temperature_steps)

	for temp_step := 0; temp_step < temperature_steps; temp_step++ {
		accept_move_func := getMoveAcceptanceFunc(temp_steps[temp_step], random)
		for step_in_temp := 0; step_in_temp < steps_per_temp; step_in_temp++ {
			current_model = current_model.GenerateMove(accept_move_func)
		}
	}
	return current_model
}

func getMoveAcceptanceFunc(temp float64, random *rand.Rand) func(log_p_curr float64, log_p_next float64) bool {
	return func(log_p_curr float64, log_p_next float64) bool {
		uniform_rand := random.Float64()
		pc := log_p_curr - math.Max(log_p_curr, log_p_next) // having some problems with numerical stability so cheating this way, TODO fix this
		pn := log_p_next - math.Max(log_p_curr, log_p_next)
		return math.Exp(-(math.Exp(pc)-math.Exp(pn))/temp) > uniform_rand
	}
}

func GetExponentialTemperatureSteps(start_temp float64, total_steps int) []float64 {
	//epsilon := 0.001
	log_start := math.Log(start_temp)
	log_end := math.Log(0.000000000001) //math.Log(epsilon)
	temps := make([]float64, total_steps)
	step := (log_end - log_start) / float64(total_steps)

	for i := 0; i < total_steps; i++ {
		temps[i] = math.Exp(log_start + (step * float64(i)))
	}
	return temps
}
