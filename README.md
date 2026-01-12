# Mathematical_Modeling_simulations

# p53-Mdm2 Stochastic Simulation Repository

This repository contains MATLAB implementations of the Gillespie Stochastic Simulation Algorithm (SSA) to model the p53 signaling network. The scripts replicate the foundational work of **Proctor & Gray (2008)** and extend it with modern regulatory motifs (miR-192 and Wip1) to study system robustness and aging.

## 📂 Repository Contents

### 1. Core ARF Model

* **Description:** Implementation of the standard p53-Mdm2-ARF negative feedback loop.
* **Features:** Simulates the system under normal homeostasis and acute DNA damage (Irradiation). Validates the necessity of the Mdm2-mRNA time delay for sustained oscillations.

### 2. Aging Simulations (Proteostasis Impairment)

* **Description:** Modifies the ARF model to simulate cellular aging.
* **Mechanism:** Simulates a ~50% decline in proteasomal efficiency by **halving the degradation rate constants** (`k_deg`) for p53, Mdm2, and ARF.

### 3. Chronic Stress Model (ROS Accumulation)

* **Description:** Simulates the gradual accumulation of Reactive Oxygen Species (ROS) over time, rather than acute irradiation.
* **Features:** Models the threshold-dependent activation of p53 and the "memory effect" caused by slow ARF turnover after repair.

### 4. Virtual Dose-Response Experiment

* **Description:** A batch simulation script designed to validate the **Digital Signaling Hypothesis**.
* **Method:** Runs 100 stochastic simulations at three distinct damage intensities (0.1, 0.5, and 5.0 Gy).
* **Analysis:** Includes signal processing to detect peaks and classify cells as "Responders" vs. "Non-Responders," confirming that cells encode damage severity in the *number* of pulses rather than amplitude.

### 5. ATM Model

* **Description:** An alternative topology where the stress signal is mediated by the ATM kinase (triggered by Double Strand Breaks) rather than ARF.
* **Features:** Includes explicit transcription steps to introduce the noise required for oscillation in this topology.

### 6. The Integrated Super-Model (ATM + miR-192 + Wip1)

* **Description:** Extending the ATM backbone with two critical regulatory loops found in recent literature (Moore et al., 2015; Luo et al., 2023):
* **Robustness Module (miR-192):** Adds a **positive feedback loop** where p53 activates miR-192.
* **Termination Module (Wip1):** Adds a coupled **negative feedback loop** where p53 activates Wip1, a phosphatase that resets ATM and p53.
