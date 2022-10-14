#pragma once

#include "../eclic.h"
#include "iss.h"
#include "nuclei_csr.h"
#include "nuclei_irq_if.h"

namespace rv32 {

struct NUCLEI_ISS : public nuclei_external_interrupt_target, public ISS {
	nuclei_csr_table csrs;
	bool clic_irq = false;
	uint32_t irq_id = 0;
	ECLIC<NUMBER_INTERRUPTS, MAX_PRIORITY>* eclic = nullptr;

	NUCLEI_ISS(uint32_t hart_id, bool use_E_base_isa = false) : ISS(hart_id, use_E_base_isa){};

	nuclei_csr_table* get_csr_table() override;
	uint32_t get_csr_value(uint32_t addr) override;
	void set_csr_value(uint32_t addr, uint32_t value) override;

	void trigger_external_interrupt(uint32_t irq_id) override;
	void clear_external_interrupt(uint32_t irq_id) override;

	void return_from_trap_handler(PrivilegeLevel return_mode) override;

	void switch_to_trap_handler();

	void run_step() override;
};

}  // namespace rv32
