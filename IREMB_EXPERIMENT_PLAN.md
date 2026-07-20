# iREMB 실험 실행 계획 (CPC 재투고 B단계)

작성일: 2026-07-20 · 상위 문서: [CPC_RESUBMISSION_PLAN.md](CPC_RESUBMISSION_PLAN.md)

**진행 순서: B(iREMB 실험) → C(원고 수정) → A(코드 공개)**
코드 공개는 원고가 어느 정도 정리된 뒤 진행한다.

---

## 0. 착수 전 해결할 불일치 ⚠️

`scripts/apply_patch.sh`가 복사하는 파일 수와 원고의 주장이 맞지 않는다.

| 항목 | `apply_patch.sh` | 원고 (`PAPER_draft.md`) |
|---|---|---|
| A-1 신규 OMP | 40개 `pair_*_omp.{h,cpp}` (= 약 20 style) | **35 style** |
| A-3 restrict backport | **15개** `pair_*.cpp` | **102개 파일** |
| A-4 transcendental | 4개 (nm/cut 계열) | 3 style (nm/cut, lj/pirani, mie/cut) |

실제 저장소: `src/OPENMP/` 79파일, `src/EXTRA-PAIR/` 6파일, `src/*.cpp` 15파일.

**CPC는 기탁 코드와 논문 기술이 일치하는지 확인하는 저널이다.** 실험 착수 전에
(a) 저장소가 논문의 전체 산출물을 담고 있는지, (b) 아니면 원고 수치가 과대인지를 확정해야 한다.
앞서 확인된 "35 vs 36종" 불일치와 함께 처리한다.

→ **이 확인이 끝나야 E2(bit-identical 검증)의 대상 목록이 정해진다.**

---

## 1. 코드 이전

저장소는 비공개지만 iREMB에 GitHub SSH 키가 등록돼 있어 clone 가능하다(HTTPS가 아닌 SSH URL 사용).

```bash
ssh iremb
cd /home/cyyou68/repos
git clone git@github.com:mirryou-maker/lammps-CO.git      # 16.6 MB (lammps-src는 .gitignore)
# LAMMPS 본체는 Ultrafast 작업용으로 이미 clone되어 있음 (동일 커밋 91d4111a94)
cp -r lammps lammps-CO-tree                                # 패치 적용용 별도 트리
bash lammps-CO/scripts/apply_patch.sh /home/cyyou68/repos/lammps-CO-tree
```

**LAMMPS 커밋이 Ultrafast-Demag와 동일한 `91d4111a94`** 이므로 재다운로드가 불필요하다.
`apply_patch.sh`가 이미 Linux용으로 존재하고, CMakeLists 수정도 불필요하다
(LAMMPS가 `RegisterStylesExt`로 `_omp.cpp`를 자동 발견).

---

## 2. 빌드 매트릭스

원고의 A-2(컴파일러 플래그) 주장을 GCC로 교차검증하려면 플래그 변형이 필요하다.
Haswell은 AVX2+FMA를 지원하나 **AVX-512는 없다**(INTEL 패키지 SVML 경로는 여전히 꺼짐).

| # | 이름 | 트리 | 플래그 | 용도 |
|---|---|---|---|---|
| 1 | `base-O3` | upstream | `-O3` | 기준선 |
| 2 | `base-avx2` | upstream | `-O3 -march=core-avx2` | A-2 기준선 |
| 3 | `co-O3` | patched | `-O3` | A-1/A-3/A-4 효과 |
| 4 | `co-avx2` | patched | `-O3 -march=core-avx2` | 주 측정 빌드 |
| 5 | `co-avx2-fast` | patched | `+ -ffast-math` | **MSVC `/fp:fast` 6.3% 주장 교차검증** |
| 6 | `co-avx2-lto` | patched | `+ -flto` | LTO 재현성 확인 (원고에서 "미미/변동") |

공통 CMake 옵션:

```
-D CMAKE_BUILD_TYPE=Release -D BUILD_SHARED_LIBS=no
-D BUILD_MPI=on -D BUILD_OMP=on
-D PKG_OPENMP=on -D PKG_KSPACE=on -D PKG_EXTRA-PAIR=on -D PKG_MANYBODY=on
-D FFT=FFTW3
```

**툴체인**: `module load CMAKE/3.31.9 DEVTOOLSET/11 OPENMPI/4.1.4.GCC8.5` +
`export OMPI_CXX=$(command -v g++) OMPI_CC=$(command -v gcc)`
(OpenMPI가 GCC 8.5로 빌드돼 있어 래퍼 컴파일러를 DEVTOOLSET/11로 재지정. FerroX에서 검증된 패턴)

빌드는 로그인 노드가 아닌 **PBS 잡(workq, ncpus=20)** 으로 제출한다.
6개 빌드 × `-j20` ≈ 1시간.

---

## 3. 실험 매트릭스

전 실험이 CPU 전용이므로 큐는 `workq`, 자원은 `select=2:ncpus=20:mpiprocs=20` (할당 40코어 준수).

### E1. MPI × OpenMP 하이브리드 스케일링 ★ 최우선

**원고 최대 공백을 메우는 신규 데이터.** 현재 전 측정이 `BUILD_MPI=OFF` 단일 프로세스다.

| 축 | 값 |
|---|---|
| 계 | 32,000 원자 FCC Ar, `lj/cut`, NVE, 500 step (원고와 동일) |
| Pure MPI | 1, 2, 4, 8, 16, 20, 40 랭크 × 1 스레드 |
| Pure OMP | 1 랭크 × 1, 2, 4, 8, 16, 20 스레드 (1노드) |
| 하이브리드 | (2×10), (4×5), (2×20), (4×10), (8×5) — 랭크×스레드 |
| 빌드 | `base-avx2`, `co-avx2` |
| 반복 | N=5, min-of-N |

추가로 **약한 스케일링**: 랭크당 원자 수 고정(32k/rank)으로 1→40랭크.

**논문값이 가장 큰 파생 결과**: 랭크가 늘면 랭크당 원자 수가 줄어 pair 분율이 떨어진다.
따라서 **OMP 최적화 이득이 MPI 랭크 수에 따라 어떻게 변하는지**를 정량화할 수 있고,
이는 "실제 HPC 환경에서 이 최적화가 유효한가"라는 리뷰어 질문에 대한 직접적 답이 된다.

`--prefix /home/APP/Software/openmpi/4.1.4` 필수 (노드 간 실행 시 `orted` 미발견 방지).

### E2. Linux/GCC bit-identical 재검증

"MSVC 특유의 아티팩트 아니냐"는 반론을 봉쇄한다.

- 대상: §0에서 확정된 전체 파일 목록 (OMP 신규 + A-3 + A-4)
- 프로토콜: 50 step NVE, `thermo_modify format float %20.15g`, step/temp/epair/etotal/press 전수 비교
- **원고보다 강화**: 스레드 수 1/2/4/8/20 **및 MPI 랭크 1/2/4/8** 조합에서 각각 확인
  → 리듀스 순서 의존성 검증. 원고는 이 축을 명시하지 않았다.
- 판정: base vs CO 마지막 자리까지 일치

### E3. 스레드 스케일링 곡선 (1→20)

원고는 1T/8T **두 점**뿐이고, 측정 환경이 하이브리드 P/E 코어 + WSL2(자인 노이즈 ±10%)였다.

- 균질 20코어 Xeon + 네이티브 Linux + 배치 스케줄러 → 노이즈 문제 해소
- 1, 2, 4, 6, 8, 10, 12, 16, 20 스레드 곡선
- Amdahl 포화 지점을 정량화 → 현재의 정성적 서술("sub-linear scaling") 대체

### E4. 장시간 NVE 에너지 드리프트 회귀

현재 검증은 **2종 × 250 step**뿐이다. 소프트웨어 저널 기준으로 부족하다.

- 전 스타일 × **10,000 step 이상** NVE, 상대 드리프트 |ΔE/E| 보고
- serial vs OMP-N 궤적 점별 최대 차이
- **CI 스크립트 형태로 배포** → "자동 수치 검증이 방법론의 성립 조건"이라는 결론을 실증

### E5. 컴파일러 플래그 교차검증 (A-2)

빌드 1~6을 32k/500step 기준으로 비교. MSVC `/fp:fast` + AVX2의 6.3% 개선이
GCC `-march=core-avx2 -ffast-math`에서도 재현되는지 확인.

### E6. 시스템 크기 의존성

원고 Limitation이 "medium-sized 4k–16k atom"이라 자인한다. 2노드면 훨씬 키울 수 있다.

- 4k / 32k / 256k / 1M 원자 × (최적 하이브리드 구성)
- Discussion "Scalability and Problem-Size Dependence" 절을 실증으로 채움

---

## 4. 소요 시간 추정

32k/500step 1런이 8스레드에서 약 4.4초(원고 실측)이므로 개별 측정은 초 단위다.

| 단계 | 런 수 | 예상 |
|---|---|---|
| 빌드 6종 | — | ~1 h |
| E1 스케일링 | ~180 | ~40 min |
| E2 bit-identical | ~500 | ~1 h |
| E3 스레드 곡선 | ~90 | ~20 min |
| E4 장시간 NVE | ~35 × 10k step | ~1 h |
| E5 플래그 | ~60 | ~15 min |
| E6 크기 스윕 | ~40 | ~40 min |
| **합계** | | **약 5~6시간** |

**하루 안쪽에 끝난다.** Ultrafast 앙상블(30시간)이나 FerroX보다 훨씬 가볍다.

---

## 5. 스케줄링 제약

- **할당은 2노드 40코어.** E1의 40랭크 측정은 2노드를 온전히 점유하므로
  Ultrafast 앙상블 잡과 **동시 실행 불가**.
- 현재 Ultrafast 잡 어레이(`243383[]`)의 서브잡 0,1이 실행 중(태스크 0~999).
  조성비 스윕 완료 후 자연스러운 공백이 생기므로 **그때 CO 작업을 투입**한다.
- workq는 잡 어레이를 무제한 동시 실행시키므로, 어레이 사용 시
  `qalter -W max_run_subjobs=2`로 제한할 것.

---

## 6. 데이터 → 원고 매핑

| 실험 | 채우는 곳 |
|---|---|
| E1 | **신규 절** "MPI–OpenMP hybrid scaling" + Limitation (i) 해소 |
| E2 | Verification 절 강화 (플랫폼 독립성) |
| E3 | Fig. A-1 스레드 스케일링 교체, "sub-linear scaling" 정량화 |
| E4 | Verification 절 + 배포 CI 스크립트 (CPC 문서화 요건) |
| E5 | A-2 절 교차검증 (MSVC → GCC) |
| E6 | Discussion "Scalability and Problem-Size Dependence" |

E2·E4 산출물은 **CPC 기탁 패키지의 일부**가 된다(재현 가능한 검증 하네스).

---

## 7. 실행 체크리스트

- [ ] §0 파일 수 불일치 확정 (35 vs 36종, 102 vs 15파일)
- [ ] iREMB로 저장소 clone (SSH URL) + `apply_patch.sh` 적용
- [ ] 빌드 6종 (PBS 잡, workq)
- [ ] E2 소규모 선행 — 1개 스타일로 bit-identical 파이프라인 검증
- [ ] E1 → E3 → E5 → E6 → E4 순 실행 (가벼운 것부터, E4는 장시간)
- [ ] 결과 CSV를 `tools/bench/results/iremb/`에 저장, 작도 스크립트 작성
