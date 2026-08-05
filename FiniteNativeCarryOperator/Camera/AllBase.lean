import Mathlib.Data.Finset.Interval
import Mathlib.Tactic

namespace FiniteNativeCarryOperator.Camera

/-- C2 is the canonical aligned chart; natural cameras start at width three. -/
def IsSupported (camera : ℕ) : Prop := 2 ≤ camera

/-- Natural half range.  C2 is handled separately by the operator. -/
def halfRange (camera : ℕ) : ℕ := camera / 2

/-- Positive radii used by a natural saturated camera. -/
def radiusSet (camera : ℕ) : Finset ℕ := Finset.Icc 1 (halfRange camera)

/-- The aligned center numbered from zero. -/
def alignedCenter (camera index : ℕ) : ℕ :=
  if camera = 2 then 4 * (index + 1) else camera * (index + 1)

/-- C2 has one seed; a natural camera has one seed for every radius. -/
def seedCount (camera : ℕ) : ℕ :=
  if camera = 2 then 1 else halfRange camera

/-- Number of centered brackets in a cutoff of `M` centers. -/
def bracketCount (camera cutoff : ℕ) : ℕ :=
  if camera = 2 then cutoff else cutoff * halfRange camera

/-- Total number of emitted coordinates used in the score normalization. -/
def coordinateCount (camera cutoff : ℕ) : ℕ :=
  seedCount camera + bracketCount camera cutoff

@[simp] theorem halfRange_two : halfRange 2 = 1 := by
  native_decide

@[simp] theorem halfRange_four : halfRange 4 = 2 := by
  native_decide

@[simp] theorem radiusSet_two : radiusSet 2 = {1} := by
  native_decide

@[simp] theorem radiusSet_four : radiusSet 4 = {1, 2} := by
  native_decide

@[simp] theorem alignedCenter_two (index : ℕ) :
    alignedCenter 2 index = 4 * (index + 1) := by
  simp [alignedCenter]

@[simp] theorem alignedCenter_of_ne_two
    {camera : ℕ} (h : camera ≠ 2) (index : ℕ) :
    alignedCenter camera index = camera * (index + 1) := by
  simp [alignedCenter, h]

@[simp] theorem bracketCount_two (cutoff : ℕ) :
    bracketCount 2 cutoff = cutoff := by
  simp [bracketCount]

@[simp] theorem bracketCount_four (cutoff : ℕ) :
    bracketCount 4 cutoff = 2 * cutoff := by
  simp [bracketCount, halfRange, Nat.mul_comm]

@[simp] theorem coordinateCount_two (cutoff : ℕ) :
    coordinateCount 2 cutoff = cutoff + 1 := by
  simp [coordinateCount, seedCount, bracketCount, Nat.add_comm]

/-- The even camera C4 contains the antipodal radius `2`; it is not C2. -/
theorem two_mem_radiusSet_four : 2 ∈ radiusSet 4 := by
  native_decide

end FiniteNativeCarryOperator.Camera
