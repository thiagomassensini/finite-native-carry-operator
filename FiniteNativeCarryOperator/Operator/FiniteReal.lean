import FiniteNativeCarryOperator.Camera.AllBase
import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic
import Mathlib.Algebra.BigOperators.Group.Finset.Basic
import Mathlib.Tactic

open scoped BigOperators

namespace FiniteNativeCarryOperator.Operator

noncomputable section

abbrev RealPlane := ℝ × ℝ

/-- Native amplitude fixed upstream by quadratic carry mass. -/
def nativeAmplitude (n : ℕ) : ℝ :=
  if n = 0 then 0 else (n : ℝ) ^ (-(1 : ℝ) / 2)

/-- Real rotating sample.  There is no free radial exponent in this state. -/
def nativeState (time : ℝ) (n : ℕ) : RealPlane :=
  let amplitude := nativeAmplitude n
  let angle := -time * Real.log n
  (amplitude * Real.cos angle, amplitude * Real.sin angle)

/-- Centered second difference at radius `r`. -/
def centeredBracket (time : ℝ) (center radius : ℕ) : RealPlane :=
  nativeState time (center - radius) -
    (2 : ℕ) • nativeState time center +
      nativeState time (center + radius)

/-- Seeds of the C2 chart or of a natural saturated camera. -/
def seedSum (camera : ℕ) (time : ℝ) : RealPlane :=
  if camera = 2 then nativeState time 1
  else ∑ radius ∈ Camera.radiusSet camera, nativeState time radius

/-- Radius sum at one aligned center. -/
def centerBracketSum (camera index : ℕ) (time : ℝ) : RealPlane :=
  if camera = 2 then
    centeredBracket time (Camera.alignedCenter 2 index) 1
  else
    ∑ radius ∈ Camera.radiusSet camera,
      centeredBracket time (Camera.alignedCenter camera index) radius

/--
Finite native operator in `R²`.  `cutoff` means exactly that many aligned
centers.  No post-bracket calibration is inserted.
-/
def finiteNativeOperator
    (camera cutoff : ℕ) (time : ℝ) : RealPlane :=
  seedSum camera time +
    ∑ index ∈ Finset.range cutoff, centerBracketSum camera index time

@[simp] theorem nativeAmplitude_zero : nativeAmplitude 0 = 0 := by
  simp [nativeAmplitude]

@[simp] theorem nativeState_zero (time : ℝ) : nativeState time 0 = 0 := by
  simp [nativeState, nativeAmplitude]

@[simp] theorem finiteNativeOperator_zero_cutoff
    (camera : ℕ) (time : ℝ) :
    finiteNativeOperator camera 0 time = seedSum camera time := by
  simp [finiteNativeOperator]

@[simp] theorem centerBracketSum_two (index : ℕ) (time : ℝ) :
    centerBracketSum 2 index time = centeredBracket time (4 * (index + 1)) 1 := by
  simp [centerBracketSum]

/-- The implementation is literally seed plus the finite sum of center channels. -/
theorem finiteNativeOperator_eq_seed_add_centers
    (camera cutoff : ℕ) (time : ℝ) :
    finiteNativeOperator camera cutoff time =
      seedSum camera time +
        ∑ index ∈ Finset.range cutoff, centerBracketSum camera index time := by
  rfl

end
end FiniteNativeCarryOperator.Operator
