import FiniteNativeCarryOperator.Operator.Score
import Mathlib.Analysis.Calculus.Deriv.Basic
import Mathlib.Tactic

namespace FiniteNativeCarryOperator.Certification

noncomputable section

/-- Closed real interval with an explicit order witness. -/
structure ClosedInterval where
  lower : ℝ
  upper : ℝ
  ordered : lower ≤ upper

namespace ClosedInterval

def Contains (I : ClosedInterval) (x : ℝ) : Prop := I.lower ≤ x ∧ x ≤ I.upper

def Width (I : ClosedInterval) : ℝ := I.upper - I.lower

theorem width_nonneg (I : ClosedInterval) : 0 ≤ I.Width := by
  exact sub_nonneg.mpr I.ordered

end ClosedInterval

/-- A sound enclosure of a scalar function on an input interval. -/
structure ScalarEnclosure (f : ℝ → ℝ) (domain : ClosedInterval) where
  range : ClosedInterval
  sound : ∀ x, domain.Contains x → range.Contains (f x)

/-- A sound coordinatewise enclosure of an `R²` function. -/
structure PlaneEnclosure
    (f : ℝ → Operator.RealPlane) (domain : ClosedInterval) where
  xRange : ClosedInterval
  yRange : ClosedInterval
  soundX : ∀ x, domain.Contains x → xRange.Contains (f x).1
  soundY : ∀ x, domain.Contains x → yRange.Contains (f x).2

/-- Proof obligation for a certified finite zero. -/
structure FiniteZeroCertificate (camera cutoff : ℕ) where
  interval : ClosedInterval
  witness : ℝ
  witness_mem : interval.Contains witness
  zero : Operator.IsFiniteNativeCarryOperatorZero camera cutoff witness
  unique : ∀ t, interval.Contains t →
    Operator.IsFiniteNativeCarryOperatorZero camera cutoff t → t = witness

/--
Proof obligation for a certified nondegenerate minimum of the finite resultant
energy.  This is distinct from a vector zero certificate.
-/
structure FiniteMinimumCertificate (camera cutoff : ℕ) where
  interval : ClosedInterval
  witness : ℝ
  witness_mem : interval.Contains witness
  stationary : deriv
      (fun t => Operator.quadraticEnergy
        (Operator.finiteNativeOperator camera cutoff t)) witness = 0
  strictLocalMinimum : ∀ t, interval.Contains t → t ≠ witness →
    Operator.quadraticEnergy
        (Operator.finiteNativeOperator camera cutoff witness) <
      Operator.quadraticEnergy
        (Operator.finiteNativeOperator camera cutoff t)

/-- Uniform norm bound relating a finite operator to a proposed limit operator. -/
structure TailBound
    (finite : ℕ → ℝ → Operator.RealPlane)
    (limit : ℝ → Operator.RealPlane)
    (domain : ClosedInterval) where
  error : ℕ → ℝ
  error_nonneg : ∀ cutoff, 0 ≤ error cutoff
  sound : ∀ cutoff t, domain.Contains t →
    Operator.quadraticEnergy (limit t - finite cutoff t) ≤ error cutoff ^ 2
  tendsToZero : Filter.Tendsto error Filter.atTop (nhds 0)

/-- Tail control for a resultant together with its first two time derivatives. -/
structure TailBoundThroughSecondDerivative
    (finite finiteFirst finiteSecond : ℕ → ℝ → Operator.RealPlane)
    (limit limitFirst limitSecond : ℝ → Operator.RealPlane)
    (domain : ClosedInterval) where
  resultant : TailBound finite limit domain
  firstDerivative : TailBound finiteFirst limitFirst domain
  secondDerivative : TailBound finiteSecond limitSecond domain

/-- Proof obligation for a unique strict minimum of a limiting energy. -/
structure LimitMinimumCertificate
    (limit : ℝ → Operator.RealPlane) where
  interval : ClosedInterval
  witness : ℝ
  witness_mem : interval.Contains witness
  stationary : deriv
      (fun t => Operator.quadraticEnergy (limit t)) witness = 0
  strictLocalMinimum : ∀ t, interval.Contains t → t ≠ witness →
    Operator.quadraticEnergy (limit witness) <
      Operator.quadraticEnergy (limit t)

/--
A cutoff-uniform residual family for one fixed limiting witness.  This is the
exact infinite obligation left after producing nested finite ledgers.
-/
structure VanishingLimitResidualCertificate
    (limit : ℝ → Operator.RealPlane) where
  witness : ℝ
  error : ℕ → ℝ
  error_nonneg : ∀ cutoff, 0 ≤ error cutoff
  energy_le : ∀ cutoff,
    Operator.quadraticEnergy (limit witness) ≤ error cutoff
  tendsToZero : Filter.Tendsto error Filter.atTop (nhds 0)

/-- A vanishing family of upper bounds forces exact real-plane vanishing. -/
theorem VanishingLimitResidualCertificate.witness_zero
    {limit : ℝ → Operator.RealPlane}
    (certificate : VanishingLimitResidualCertificate limit) :
    limit certificate.witness = 0 := by
  apply (Operator.quadraticEnergy_eq_zero_iff (limit certificate.witness)).mp
  apply le_antisymm
  · by_contra hnot
    have hpositive : 0 < Operator.quadraticEnergy (limit certificate.witness) :=
      lt_of_not_ge hnot
    have heventually : ∀ᶠ cutoff in Filter.atTop,
        certificate.error cutoff <
          Operator.quadraticEnergy (limit certificate.witness) :=
      certificate.tendsToZero.eventually (Iio_mem_nhds hpositive)
    rcases heventually.exists with ⟨cutoff, hstrict⟩
    exact (not_lt_of_ge (certificate.energy_le cutoff)) hstrict
  · exact Operator.quadraticEnergy_nonneg (limit certificate.witness)

/-- The simple polynomial envelope used for the uniform C3 tail witness. -/
def polynomialTailEnvelope (constant : ℝ) (cutoff : ℕ) : ℝ :=
  constant * (1 / ((cutoff : ℝ) + 1)) ^ 5

/-- Every fixed multiple of the polynomial tail envelope tends to zero. -/
theorem polynomialTailEnvelope_tendsToZero (constant : ℝ) :
    Filter.Tendsto (polynomialTailEnvelope constant)
      Filter.atTop (nhds 0) := by
  change Filter.Tendsto
    (fun cutoff : ℕ => constant * (1 / ((cutoff : ℝ) + 1)) ^ 5)
    Filter.atTop (nhds 0)
  simpa [one_div, inv_pow] using
    (tendsto_const_nhds (x := constant)).mul
      ((tendsto_one_div_add_atTop_nhds_zero_nat (𝕜 := ℝ)).pow 5)

/-- A nonnegative error below the polynomial witness also tends to zero. -/
theorem tendsto_zero_of_le_polynomialTailEnvelope
    (constant : ℝ) (error : ℕ → ℝ)
    (error_nonneg : ∀ cutoff, 0 ≤ error cutoff)
    (error_le : ∀ cutoff, error cutoff ≤
      polynomialTailEnvelope constant cutoff) :
    Filter.Tendsto error Filter.atTop (nhds 0) := by
  exact squeeze_zero error_nonneg error_le
    (polynomialTailEnvelope_tendsToZero constant)

/--
Localization of one stationary witness from approximate stationary centers.
The slope-control field is the exact mean-value inequality supplied by a
strict positive lower bound for the limiting stationary derivative.
-/
structure StationaryLocalizationCertificate where
  stationary : ℝ → ℝ
  witness : ℝ
  stationary_witness : stationary witness = 0
  center : ℕ → ℝ
  slopeLower : ℝ
  slopeLower_pos : 0 < slopeLower
  stationaryError : ℕ → ℝ
  stationaryError_nonneg : ∀ cutoff, 0 ≤ stationaryError cutoff
  stationaryError_tendsToZero :
    Filter.Tendsto stationaryError Filter.atTop (nhds 0)
  centerError : ∀ cutoff,
    |stationary (center cutoff)| ≤ stationaryError cutoff
  slopeControl : ∀ cutoff,
    |witness - center cutoff| * slopeLower ≤
      |stationary witness - stationary (center cutoff)|

namespace StationaryLocalizationCertificate

/-- Radius obtained by dividing stationary-equation error by the slope margin. -/
def radius (certificate : StationaryLocalizationCertificate)
    (cutoff : ℕ) : ℝ :=
  certificate.stationaryError cutoff / certificate.slopeLower

/-- The limiting stationary witness lies inside every derived radius. -/
theorem witness_distance_le_radius
    (certificate : StationaryLocalizationCertificate) (cutoff : ℕ) :
    |certificate.witness - certificate.center cutoff| ≤
      certificate.radius cutoff := by
  apply (le_div_iff₀ certificate.slopeLower_pos).2
  calc
    |certificate.witness - certificate.center cutoff| *
        certificate.slopeLower ≤
      |certificate.stationary certificate.witness -
        certificate.stationary (certificate.center cutoff)| :=
      certificate.slopeControl cutoff
    _ = |certificate.stationary (certificate.center cutoff)| := by
      rw [certificate.stationary_witness]
      simp only [zero_sub, abs_neg]
    _ ≤ certificate.stationaryError cutoff :=
      certificate.centerError cutoff

/-- Vanishing stationary-equation error produces vanishing localization radii. -/
theorem radius_tendsToZero
    (certificate : StationaryLocalizationCertificate) :
    Filter.Tendsto certificate.radius Filter.atTop (nhds 0) := by
  change Filter.Tendsto
    (fun cutoff => certificate.stationaryError cutoff /
      certificate.slopeLower) Filter.atTop (nhds 0)
  simpa using
    certificate.stationaryError_tendsToZero.div_const certificate.slopeLower

/-- A fixed velocity cap times the localization radius also tends to zero. -/
theorem velocity_mul_radius_tendsToZero
    (certificate : StationaryLocalizationCertificate) (velocityCap : ℝ) :
    Filter.Tendsto (fun cutoff => velocityCap * certificate.radius cutoff)
      Filter.atTop (nhds 0) := by
  simpa using
    (tendsto_const_nhds (x := velocityCap)).mul certificate.radius_tendsToZero

end StationaryLocalizationCertificate

/--
A stationary residual split into its corrected-center, oriented-tail, and
localization contributions.  Recording the three limits separately prevents a
finite tail estimate from hiding the still-open corrected-center obligation.
-/
structure DecomposedVanishingLimitResidualCertificate
    (limit : ℝ → Operator.RealPlane) where
  witness : ℝ
  coreError : ℕ → ℝ
  tailError : ℕ → ℝ
  localizationError : ℕ → ℝ
  coreError_nonneg : ∀ cutoff, 0 ≤ coreError cutoff
  tailError_nonneg : ∀ cutoff, 0 ≤ tailError cutoff
  localizationError_nonneg : ∀ cutoff, 0 ≤ localizationError cutoff
  energy_le : ∀ cutoff,
    Operator.quadraticEnergy (limit witness) ≤
      (coreError cutoff + tailError cutoff + localizationError cutoff) ^ 2
  coreError_tendsToZero : Filter.Tendsto coreError Filter.atTop (nhds 0)
  tailError_tendsToZero : Filter.Tendsto tailError Filter.atTop (nhds 0)
  localizationError_tendsToZero :
    Filter.Tendsto localizationError Filter.atTop (nhds 0)

/-- Vanishing of all three residual components forces the limiting vector zero. -/
theorem DecomposedVanishingLimitResidualCertificate.witness_zero
    {limit : ℝ → Operator.RealPlane}
    (certificate : DecomposedVanishingLimitResidualCertificate limit) :
    limit certificate.witness = 0 := by
  apply VanishingLimitResidualCertificate.witness_zero
    { witness := certificate.witness
      error := fun cutoff =>
        (certificate.coreError cutoff + certificate.tailError cutoff +
          certificate.localizationError cutoff) ^ 2
      error_nonneg := fun cutoff => sq_nonneg _
      energy_le := certificate.energy_le
      tendsToZero := by
        simpa using
          ((certificate.coreError_tendsToZero.add
            certificate.tailError_tendsToZero).add
            certificate.localizationError_tendsToZero).pow 2 }

/-- Final bridge expected from future interval arithmetic and a rigorous tail bound. -/
structure LimitZeroCertificate
    (limit : ℝ → Operator.RealPlane) where
  interval : ClosedInterval
  witness : ℝ
  witness_mem : interval.Contains witness
  zero : limit witness = 0
  unique : ∀ t, interval.Contains t → limit t = 0 → t = witness

end
end FiniteNativeCarryOperator.Certification
