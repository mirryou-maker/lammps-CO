# 검증 결과와 발견된 결함 (2026-07-21)

iREMB에서 원본 LAMMPS와 A-3 패치본을 같은 커밋(`91d4111a9`)에서 빌드하고,
LAMMPS 자체 `unittest/force-styles` 회귀 스위트로 대조한 결과.
상위 문서: [CPC_RESUBMISSION_PLAN.md](CPC_RESUBMISSION_PLAN.md), [IREMB_EXPERIMENT_PLAN.md](IREMB_EXPERIMENT_PLAN.md)

---

## 1. 파일 수 확정

| 항목 | 원고(수정 전) | 실측 |
|---|---|---|
| A-3 수정 파일 | 102 (Table S6의 카테고리 어림 합계) | **228** |
| 그중 되돌린 파일 | — | 2 (아래 §3) |
| **최종 A-3 파일** | — | **226** |
| 신규 OMP style | 35 | **36** (+ 기존 4개 갱신 = 배포 40파일) |
| 배포 저장소가 담고 있던 A-3 파일 | **15** | 228로 동기화 |

원고의 "102"는 REPORT.md와 Table S6 모두 `~102`로 물결표를 단 **카테고리별 추정치**였는데,
본문에서 확정 수치로 승격되어 있었다. 실제 파일 목록은 `Paper/tableS6_a3_files.csv`
(파일 경로 · 패키지 · `pair_style`명 · `_noalias`/`fxtmp`/`dbl3_t` 보유 여부, 38개 패키지).

**배포본이 15개뿐이었으므로, 공개 저장소만으로는 논문의 A-3 결과를 재현할 수 없었다.**

---

## 2. 컴파일 불가 결함 10건

전부 동일한 결함: A-3 변환이 `int inum = list->inum;`을 추가하면서
기존 선언줄(`int i,j,ii,jj,inum,jnum,...;`)에서 `inum` 제거를 빠뜨림 → `redeclaration of 'int inum'`.

| 패키지 | 파일 |
|---|---|
| ASPHERE | `pair_gayberne`, `pair_resquared`, `pair_ylz` |
| COLLOID | `pair_brownian`, `pair_brownian_poly` |
| GRANULAR | `pair_gran_hertz_history`, `pair_gran_hooke`, `pair_gran_hooke_history` |
| MISC | `pair_tracker` |
| SPH | `pair_sph_heatconduction` |

**원인**: Windows/MSVC 빌드에 이 패키지들이 포함되지 않아 **한 번도 컴파일된 적이 없다.**
원고의 검증 1단계("컴파일, warnings-as-errors, 신규 경고 없음")는 이 파일들을 통과시킨 것이 아니라
아예 대상에 넣지 않았다. 생성 시점 워크플로가 이 파일들에 대해 **어떤 피드백 신호도 받지 못한 것**이며,
모델의 문제가 아니라 하네스의 문제다.

전부 수정 완료(정상 파일과 동일한 패턴으로 기존 선언줄에서 해당 이름 제거).

> 참고: 정적 스캐너는 18개 파일을 지목했으나 검증 결과 **9개는 오탐**이었다.
> `pair_soft.cpp`처럼 루프 내부에서 재선언하는 경우는 shadowing이라 합법이다.
> 휴리스틱이 아니라 컴파일러 판정(`make -k`로 전수 수집)을 근거로 삼아야 한다.

---

## 3. 수치 결함

### 3.1 대수 오류 1건 — `EXTRA-PAIR/pair_lj_expand_sphere.cpp`

```cpp
- fpair = factor_lj * forcelj * rshift / r;   // 원본
+ fpair = factor_lj * forcelj / rshift / r;   // A-3 변환 후
```

곱셈이 나눗셈으로 전사되어 힘이 `rshift²`배 틀어졌다.
회귀 스위트가 즉시 검출(**오차 3.8~5.3, 허용치 7.5×10⁻¹³**). 수정 완료.

**이것이 원고의 "bit-identical verification for every modification" 주장에 대한 직접적 반례다.**
전수 검증이 실제로 수행되었다면 반드시 걸렸어야 한다.

### 3.2 부동소수점 재결합 3건 → 되돌림

| 파일 | 검출 경로 | 오차 |
|---|---|---|
| `src/pair_buck_coul_cut.cpp` | `buck_coul_cut_qeq_point` 테스트 | 7.1×10⁻¹² (허용 3.75×10⁻¹²) |
| `src/INTERLAYER/pair_kolmogorov_crespi_z.cpp` | 〃 | 〃 |
| `src/MANYBODY/pair_polymorphic.cpp` | 전용 검증 덱 (§5) | 압력 마지막 2자리 (~1×10⁻¹²) |

O(1) 오류가 아니라 누적 순서 변경에 따른 반올림 차이지만, 허용치를 넘으므로 bit-identical이 아니다.
**세 파일 모두 원본으로 되돌렸다.** `pair_buck_coul_cut.cpp`는 원래 배포되던 15개 중 하나였으므로,
이 편차는 공개 저장소에 계속 존재해 왔다.

→ **최종 A-3 파일 수: 228 − 3 = 225**

---

## 4. upstream LAMMPS 버그 1건 — 발견 및 수정 (본 연구의 부수 기여)

`ManybodyPairStyle:kolmogorov_crespi_z`가 패치본에서만 실패했다. 실패한 하위 테스트가
`PairStyle.omp` / `PairStyle.extract_omp`여서 처음에는 신규 OMP 변형(A-1)의 결함으로 보였으나,
**최소 입력으로 재현한 결과 수정하지 않은 upstream 빌드에서도 동일하게 실패했다.**

### 원인

upstream은 `hybrid/overlay`의 OMP 별칭을 같은 클래스에 등록한다:

```cpp
// src/pair_hybrid_overlay.h
PairStyle(hybrid/overlay,PairHybridOverlay);
PairStyle(hybrid/overlay/omp,PairHybridOverlay);
```

따라서 `-sf omp`로 실행하면 `Force::new_pair()`가 `hybrid/overlay/omp`를 찾아 `sflag=1`이 되고,
`force->pair_style`이 **`"hybrid/overlay/omp"`** 로 저장된다. 그런데 자식 스타일의 검사는
정확 일치를 요구한다:

```cpp
// src/INTERLAYER/pair_kolmogorov_crespi_z.cpp:206
if (strcmp(force->pair_style, "hybrid/overlay") != 0)
    error->all(FLERR, "ERROR: requires hybrid/overlay pair_style");
```

결과적으로 **`pair_style hybrid/overlay kolmogorov/crespi/z`는 `-sf omp` 환경에서 무조건 실패한다.**
upstream 테스트가 이를 잡지 못한 이유는, OMP 변형이 없으면 해당 하위 테스트가 skip되기 때문이다.
**본 연구가 OMP 변형을 추가하면서 비로소 이 경로가 실행되어 버그가 드러났다.**

### 수정 (1줄)

```cpp
-  if (strcmp(force->pair_style, "hybrid/overlay") != 0)
+  if (!utils::strmatch(force->pair_style, "^hybrid/overlay"))
```

수정 후 **325/325 전부 통과**. 이 변경은 A-3 변환과 무관한 별개의 upstream 버그 수정이며,
LAMMPS 본체에 별도 PR로 제출할 가치가 있다.

동일한 정확 일치 패턴이 `FEP/compute_fep.cpp:222`, `FEP/fix_adapt_fep.cpp:308`에도 있어
같은 조건에서 오작동할 수 있다(미검증, 확인 권장).

---

## 5. 최종 검증 상태

원본과 패치본을 동일 커밋에서 빌드하고 `most` 패키지 프리셋(외부 라이브러리 3종은
로그인 노드에서 사전 다운로드 — 계산 노드는 네트워크 차단)으로 325개 pair-style 테스트 대조:

| 빌드 | 결과 |
|---|---|
| base (원본) | **325/325 통과** |
| patched (A-3 226 + A-1 40) | **324/325 통과** — 잔여 1건은 §4 |

초기 축소 패키지 구성(262 테스트)에서는 base·patched 결과가 **완전히 일치(차이 0건)** 했다.

### 검증 커버리지 — 완결

`most` 패키지 구성으로 재측정한 결과, 공백은 46개가 아니라 **11개**였다
(초기 46은 축소 패키지 빌드(262 테스트) 기준이었고, `examples/` 스캔이 `.lmp` 확장자와
hybrid 하위 스타일을 놓치고 있었다).

| 경로 | style 수 |
|---|---|
| force-style 회귀 스위트 (325 테스트) | 대다수 |
| `examples/` 기존 입력덱 | 2 (`sph/idealgas`, `sph/taitwater/morris`) |
| **전용 검증 덱 신규 작성** | **11** |
| 미검증 | **0** |

신규 작성한 11개는 `tools/verify/`에 하네스와 함께 배포한다.

| style | 셋업 |
|---|---|
| `brownian`, `brownian/poly` | `atom_style sphere`; `brownian/poly`는 `newton off` 필요 |
| `dsmc` | |
| `eam/apip` | `atom_style apip`; **표준 EAM setfl 파일(`Cu_smf7.eam`)을 그대로 수용** |
| `multi/lucy/rx`, `table/rx` | `examples/PACKAGES/dpd-react/dpdrx-shardlow`의 `fix rx` 셋업을 재활용. `table/rx`용 `rxn.table`은 생성 |
| `polymorphic` | `potentials/CuTa_eam.poly` |
| `sph/heatconduction`, `sph/lj` | `atom_style sph` |

`eam/apip`은 배포 퍼텐셜 파일이 없어 검증 불가로 예상했으나, 표준 EAM 파일로 동작해
**미검증 style이 하나도 남지 않았다.**

비교 방식: 동일 커밋에서 빌드한 원본·패치본으로 같은 덱을 실행하고 thermo 출력을 문자 단위로 대조.
벽시계 시간·구간별 타이밍 표·호스트 배너만 제외하고 나머지는 전부 일치해야 한다.

---

## 6. 원고 반영 현황

| 위치 | 상태 |
|---|---|
| §A-3 본문 — 102 → 228, 제외 기준 삭제, Table S6 참조 | ✅ 반영 |
| Correctness Verification — 2차 캠페인·결함 10건 기술 | ✅ 반영 |
| Error Analysis — "98 of 102 (96%)" → "218 of 228 (95.6%)" | ✅ 반영 |
| Abstract·Figure 1 캡션·Supplementary Table S2 설명의 102 | ⬜ 미반영 |
| §4의 A-1 결함 반영 | ⬜ 미해결 사항 확정 후 |
| Table S6 본체를 CSV 목록으로 교체 | ⬜ |

---

## 7. 방법론적 함의 (원고에 쓸 가치가 있는 내용)

LLM이 생성한 228개 파일에서 **컴파일 오류 10건(4.4%)과 수치 오류 1건**이 나왔고,
자동 회귀 검증이 그 전부를 잡아냈다. 이는 원고의 결론 —
"전문가 검토 없이 잘못된 생성물을 기각할 자동 수치 검증이 방법론의 성립 조건" —
을 **약화시키는 것이 아니라 실증하는 사례**다.

동시에 두 가지 하네스 요구사항을 드러낸다.

1. **생성 대상 파일이 속한 패키지 전체를 컴파일해야 한다.** 벤치마크에 쓰는 패키지만
   빌드하면, 건드린 파일 중 일부는 아무 피드백 없이 통과된 것처럼 집계된다.
2. **파일 수를 세어야 한다.** 카테고리별 어림값으로는 "건드렸지만 검증되지 않은 파일"이
   생겨도 드러나지 않는다.
