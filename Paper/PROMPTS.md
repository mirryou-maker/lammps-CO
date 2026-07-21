# LAMMPS 최적화 프로젝트 — 사용자 프롬프트 모음

> **프로젝트**: LLM-Assisted Systematic Optimization of the LAMMPS MD Simulator  
> **기간**: 2026-06-07 ~ 2026-06-23  
> **총 세션**: 2개 (Session-1: 2026-06-07, Session-2: 2026-06-11~23)  
> **총 프롬프트**: 47개 (중복 제거)

---

## Phase 0 — 프로젝트 시작 및 분석 계획

### Session-1 (2026-06-07)

**[P-01]** `2026-06-07`  
```
이 프로젝트는 Molecular dynamics 시뮬레이션 프로그램인 LAMMPS (https://www.lammps.org)의
코드를 분석하고, Claude code를 이용한 최적화를 통해서 코드의 성능을 향상시키고자 해.

먼저 다음 사항을 순차적으로 확인해주고 그 결과를 Claude.md 에 정리해서 저장해줘.
1. https://www.lammps.org 에서 소스코드와 manual, document 들을 읽어서 LAMMPS의 구조를 분석해줘.
2. 소스코드를 분석해서 Claude-code를 이용해서 개선의 가능성이 있다면 개선할 수 있는 항목들을
   검토해서 개선 계획을 세워줘.
3. 세워진 계획들을 구현할 경우, 어느 정도의 속도 개선이 가능할지 예측해줘.
4. 분석 결과와 개선 계획을 Claude.md에 저장해줘.
```
> → CLAUDE.md 작성 (LAMMPS 구조 분석, A-1~A-3/B-1~B-2 개선 계획, 속도 예측)

**[P-02]** `2026-06-07`  
```
현재 분석 내용을 기억해두고, 다시 개발을 시작할 때 권장작업을 추천해줘.
```
> → 메모리 저장 및 재시작 권고사항 정리

---

## Phase 1 — 구현 전략 수립 및 B-1 벤치마크 하네스

### Session-2 (2026-06-11)

**[P-03]** `2026-06-11`  
```
LAMMPS를 최적화하기 위한 전략을 step by step으로 제시해줘.
```
> → 5단계 구현 계획 (B-1 벤치마크 → B-2 가속 진단 → A-2 빌드 플래그 → A-1 OMP 포팅 → A-3 백포트)

**[P-04]** `2026-06-11`  
```
Step1을 구현해줘.
```
> → `tools/bench/Run-Benchmark.ps1`, `LammpsLog.psm1`, `Build-Lammps.ps1` 구현; 기준선 측정 (Loop time ≈ 0.29s)

**[P-05]** `2026-06-11`  
```
여기까지 저장하고, 다음 시작시 계속 진행할꺼야.
```
> → 메모리 저장

---

## Phase 2 — B-2 가속 패키지 진단 + Step 3 빌드 플래그

**[P-06]**  
```
먼저 step 2(가속 패키지..)를 먼저 처리 후 Step3의 진행을 물어봐줘.
```
> → `tools/accel/Check-AccelCoverage.ps1` 구현; OPENMP 빌드 검증 (4t: 3.51×, 8t: 4.18×)

**[P-07]**  
```
step1을 먼저 실행 후 step3의 실행여부를 물어봐줘.
```
> → 기준선 재측정 (재현성 확인)

**[P-08]** `2026-06-14`  
```
step 3를 순차적으로 구현해줘.
```
> → `/arch:AVX2` (+8.4%), `/fp:fast` (+3.1%), LTO (−3.6%), OMP+AVX2 비교 벤치마크 완료

**[P-09]** `2026-06-14`  
```
다음 추천단계는?
```
> → A-1 (OMP 변형 포팅) 권고

---

## Phase 3 — A-1: OpenMP pair_style 변형 포팅 (36개)

**[P-10]**  
```
step 4를 수행후 결과를 보고하고 Step 5의 실행여부를 물어봐줘.
```
> → 첫 OMP 포팅: `pair_lj_class2_soft_omp` (bit-identical, 4t: 3.12×)

**[P-11]**  
```
중단된 작업을 계속 실행해줘
```
> → `pair_nm_cut_split_omp` 포팅 (4t: 2.97×)

**[P-12]**  
```
하던 작업을 계속 진행해줘.
```
> → `pair_morse_soft_omp` 포팅 (4t: 3.18×)

**[P-13]**  
```
현재 상태를 보고하고, 다음 작업을 추천해줘.
```
> → 포팅 진행상황 보고 (3개 완료)

**[P-14]** `2026-06-17`  
```
다시 시작해줘.
```
> → 포팅 재개

**[P-15]** `2026-06-17`  
```
추천한 ROI가 높은 순서대로 순차적으로 진행해줘.
```
> → `coul/cut/soft/gapsys`, `born/coul/dsf`, `born/coul/dsf/cs`, `born/coul/wolf/cs` 등 대량 포팅 진행

**[P-16]** `2026-06-18`  
```
남은 추천 사항들은?
```
> → 남은 포팅 목록 (36개 중 ~20개 완료 시점)

**[P-17]** `2026-06-18`  
```
B-1을 먼저 진행 후 A-2까지 수행해주고 다음 작업을 추천해줘.
```
> → 포팅 계속 + 32k-atom 정식 벤치마크 (4t: 4.31×, 8t: 4.93×)

**[P-18]** `2026-06-18`  
```
b2 진행후 다음 작업을 추천해줘
```
> → 포팅 계속 (lj/charmmfsw, lj/class2/soft 계열, coul/* 계열)

**[P-19]** `2026-06-18`  
```
1순위 작업을 수행 후 다음 작업을 추천해줘
```
> → 포팅 계속

**[P-20]** `2026-06-18`  
```
다음 작업은?
```
> → 포팅 마무리 단계 (36개 완료)

**[P-21]** `2026-06-18`  
```
현재까지 결과를 요약해주고 다음 추천 작업을 알려줘
```
> → A-1 완료 보고 (36개, 51%→63% 커버리지), A-3 시작 권고

---

## Phase 4 — A-3: `__restrict__` / fxtmp 백포트 (15개 파일)

**[P-22]** `2026-06-18`  
```
step5를 마저 작업해줘
```
> → `pair_lj_cut`, `pair_lj_cut_coul_cut`, `pair_morse` 등 15개 파일에 `__restrict__`/`_noalias`/`fxtmp` 패턴 적용

**[P-23]** `2026-06-18`  
```
다시 전체적으로 최적화 할수 있는 부분이 있는지 확인하고, 없으면 현재까지 개선된 내용과
벤치마킹 결과를 보고서로 작성해줘
```
> → 전체 307개 `pair_*.cpp` 재스캔 (236개 최적화, 37개 Skip — ML/많은몸 포텐셜 등), REPORT.md 작성

---

## Phase 5 — 논문 작성 (npj Computational Materials)

**[P-24]** `2026-06-18`  
```
다음과 같은 내용으로 npj Computational Materials에 투고할 논문을 작성하려고해.
아래 내용을 검토해주고 수정/보완 의견을 제시해줘.
(그림은 5, 6개 이상 그려주고, 참고문헌은 30개 이상을 링크나 DOI와 함꼐 제시해줘)

초록
1. MD 시뮬레이션이 물질 연구에 중요한 점과 LAMMPS의 역할, 중요성 등을 설명
2. Claude Code (LLM)를 이용한 코드 최적화 내용 설명
...
```
> → `PAPER_draft.md` 초안 작성 (Abstract, Introduction, Results, Discussion, Methods), 52개 참고문헌, Fig.1~6 스크립트 작성

**[P-25]** `2026-06-18`  
```
위의 모든 제안들을 반영하는 것으로 계획을 세워주고, npj CM에 제출하기 위해서
제시한 조건을 만족하기 위한 아이디어가 있으면 제시해줘.
```
> → 6개 보강 항목 확정 (FEP 사례 연구, 32k 벤치마크, 통계적 유의성 등)

**[P-26]** `2026-06-18`  
```
현재까지의 내용을 기억하고, 다음 시작시 다음 작업을 추천해줘.
```
> → 메모리 저장

**[P-27]** `2026-06-18`  
```
Phase 1 부터 시작하고, 완료 후 후속작업을 제시해줘.
```
> → 32k-atom 정밀 벤치마크, 통계적 유의성 검증, Fig.1~6 재생성 (13pt 폰트, 300 DPI)

**[P-28]** `2026-06-18`  
```
Phase 2에서 option A 로 진행해줘.
```
> → CG FEP (Coarse-Grained Li⁺ desolvation) 40 ps/window 실행 결정

---

## Phase 6 — FEP 사례 연구 (OPLS-AA, 200 ps/window)

**[P-29]**  
```
A. 단계 (문서/코드)를 먼저 실행해주고, B. 재료 시뮬레이션의 option A를 위한
구체적 전략과 순서를 제시해줘.
```
> → OPLS-AA 전원자 FEP 파이프라인 설계 (EC/PC/DME, 63 λ-윈도우)

**[P-30]**  
```
우선순위 1인 python 환경 설정후에 그림 생성을 진행하고 우선순위 2를 실행하기 위한
step by step 작업을 제시해줘.
```
> → `analyze_fep_opls.py` 작성, TI+BAR 분석

**[P-31]** `2026-06-18`  
```
3번 선택 사항인 OPLS-AA 전원자 모델로 FEP 확장하고 Fig.7을 업그레이드 해주고,
supplementary를 작성해줘.
```
> → OPLS-AA 40ps FEP 실행, Fig.7 생성, `Supplementary_Materials.md` 작성

**[P-32]** `2026-06-18`  
```
현재까지 내용을 저장해주고, 시뮬레이션이 끝나면 필요한 업데이트를 수행해주고,
향후 작업을 추천해줘.
```
> → 메모리 저장, 200 ps/window 연장 실행 예약

**[P-33]** `2026-06-18`  
```
논문 품질 향상을 위한 6. 정량적 검증 강화를 위한 참고문헌과 값을 찾아서 비교해주고,
Production Run 연장을 수행해줘.
```
> → 문헌값 비교 (Wan 2016 −91.3, Takeuchi 2009 −111.8 kcal/mol), 200 ps 연장 실행 시작

**[P-34]** `2026-06-18`  
```
현재 상태는?
```
> → 200 ps FEP 진행 상황 보고 (EC 완료: −183.7±0.1, PC/DME 진행 중)

---

## Phase 7 — 논문 최종화 및 GitHub 배포

**[P-35]** `2026-06-19`  
```
npj Computational Materials에 투고하기 위한 cover letter를 작성해주고
(이 논문의 중요성과 npj Computational Materials에 적합한 논문임을 설득력있게 강조)
현재 그림들의 font가 작고 겹치는 등 가독성이 많이 떨어지는데,
그림들을 가독성을 높여서 다시 그려줘.
```
> → `Cover_Letter.md` 작성, Fig.1~7 전면 재작성 (13pt 폰트, 300 DPI)

**[P-36]** `2026-06-19`  
```
여기까지 작업을 기억해줘.
```
> → 세션 메모리 저장 (FEP 최종값 EC/PC/DME, 39개 참고문헌 정리)

**[P-37]** `2026-06-19`  
```
Ref 19/43 중복(둘 다 Chen 2021 Codex) 해결 / 미인용 참고문헌(29-37, 40, 42, 45) 처리 /
를 진행하고 다시 메모리를 업데이트해줘.
```
> → 참고문헌 52개 → 39개로 정리, `tools/fix_refs.py` 작성, 인용 번호 전면 재매핑

**[P-38]** `2026-06-19`  
```
다음 작업을 수행해줘.
1. repository URL을 github를 이용해서 생성하고 paper_draft.md의 해당 부분을 수정해줘.
2. 저자이름은 Chun-Yeol You, Affiliation은 Department of Physics and Chemistry,
   DGIST (Daegu Gyeongbuk Institute of Science and Technology)...
```
> → GitHub repo 생성 (`mirryou-maker/lammps-CO`), `tools/` 업로드, `README.md` 작성

**[P-39]** `2026-06-19`  
```
github에서 수정된 source code를 사용자가 다운 받을 수 있도록 해주고,
일단은 private 으로 설정해줘.
```
> → repo private 설정, v1.0.0 Release 생성

**[P-40]** `2026-06-19`  
```
참고문헌 38번의 저자이름이 잘못 된 것 같아.
DOI를 확인해서 저자이름을 수정하고, 본문의 표의 저자이름도 수정해줘.
```
> → CrossRef API로 DOI 10.1021/jp3011487 확인 → "Takeuchi, M. et al." 로 수정

**[P-41]** `2026-06-19`  
```
github에서 다운 받은 소스에는 최적화로 수정된 source code (*.cpp)이 하나도 없는데?
```
> → 72개 OMP 파일 + 15개 A-3 파일 → GitHub `src/OPENMP/`, `src/` 에 추가 업로드

**[P-42]** `2026-06-19`  
```
사용자가 우리가 최적화한 LAMMPS code를 다운 받아서 컴파일해서 사용할 수 있도록
github에 필요한 파일들을 모두 올려주고 설명서를 작성해줘.
```
> → `INSTALL.md` (6단계 설치 가이드), `scripts/apply_patch.sh`, `scripts/apply_patch.ps1` 작성 후 push

**[P-43]** `2026-06-19`  
```
현재까지 결과를 저장해줘.
```
> → 메모리 전면 업데이트 (GitHub 배포 완료 상태 반영)

---

## Phase 8 — PPT 제작

**[P-44]** `2026-06-22`  
```
현재까지 작성된 논문 내용을 기반으로 PPT를 만들고 싶어. 방법을 안내해줘.
```
> → Marp / PowerShell COM / python-pptx 옵션 안내

**[P-45]** `2026-06-22`  
```
python-pptx 를 설치하고 위 내용으로 PPT를 만들고 싶어
```
> → Anaconda Python 확인, `python-pptx` 설치, `tools/plots/make_ppt.py` 작성 (11슬라이드)

**[P-46]** `2026-06-22`  
```
생성된 PPT에 LAMMPS 관련 그림과 논문을 위해서 생성한 그림들 추가해줘.
```
> → Fig.1~7 + LAMMPS domain-decomp 다이어그램 삽입, 16슬라이드로 재구성

**[P-47]** `2026-06-23`  
```
내가 LAMMPS의 최적화를 입력한 prompt들을 정리해서 *.md 파일로 저장해줘.
```
> → 이 파일 (`PROMPTS.md`) 생성

---

## 프롬프트 패턴 분석

| 유형 | 개수 | 예시 |
|------|------|------|
| 단계 실행 지시 | 16 | "Step1을 구현해줘", "step 4를 수행해줘" |
| 상태 확인 / 진행 재개 | 10 | "현재 상태는?", "하던 작업을 계속 진행해줘" |
| 메모리 저장 | 5 | "여기까지 저장하고", "기억해줘" |
| 콘텐츠 생성 | 8 | "논문 작성", "PPT 만들기", "cover letter" |
| 수정 / 검증 요청 | 5 | "저자이름 수정", "github 파일 추가" |
| 전략 / 추천 질문 | 3 | "추천 작업을 알려줘", "ROI 높은 순서로" |

### 핵심 패턴

- **짧고 직관적**: 대부분 1~2문장 ("Step1을 구현해줘", "다음 작업은?")
- **순차적 진행**: "완료 후 후속 작업 추천" 패턴이 반복됨
- **맥락 위임**: 상세 조건을 CLAUDE.md/메모리에 위탁하고 프롬프트는 간결하게 유지
- **중간 체크포인트**: 긴 작업마다 "저장해줘" / "기억해줘" 로 상태 보존

---

## 산출물 요약

| 산출물 | 파일/위치 |
|--------|-----------|
| LAMMPS 분석 + 계획 | `CLAUDE.md` |
| 벤치마크 하네스 | `tools/bench/Run-Benchmark.ps1` 등 |
| 가속 진단 도구 | `tools/accel/Check-AccelCoverage.ps1` |
| A-1: 36개 OMP 변형 | `src/OPENMP/pair_*_omp.{h,cpp}` |
| A-3: 15개 restrict 최적화 | `src/pair_*.cpp` |
| FEP 파이프라인 | `tools/fep/` |
| 논문 초안 | `PAPER_draft.md` (39 refs, 7 figs) |
| 그림 7개 | `tools/plots/figure{1-7}_*.{pdf,png}` |
| Cover letter | `Cover_Letter.md` |
| 설치 가이드 | `INSTALL.md` |
| 패치 스크립트 | `scripts/apply_patch.{sh,ps1}` |
| PPT (16슬라이드) | `tools/plots/lammps_llm_optimization.pptx` |
| 이 파일 | `PROMPTS.md` |
