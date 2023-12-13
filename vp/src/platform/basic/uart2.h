#ifndef RISCV_ISA_UART2_H
#define RISCV_ISA_UART2_H

#include <tlm_utils/simple_target_socket.h>

#include <systemc>
#include <cstdlib>
#include <cstring>

#include "core/common/irq_if.h"

#define FIFO_FAIL 0
#define FIFO_SUCCESS 1
#define FIFO_SIZE 8

struct FIFO {
	uint8_t data[FIFO_SIZE];
	uint8_t read;
	uint8_t write;
	uint8_t size;

	uint8_t push(uint8_t byte) {
		if (size == FIFO_SIZE)
			return FIFO_FAIL;

		data[write++] = byte;
		write %= FIFO_SIZE;
		size++;
		return FIFO_SUCCESS;
	}

	uint8_t pop(uint8_t *pByte) {
		if (size == 0)
			return FIFO_FAIL;

		*pByte = data[read++];
		read %= FIFO_SIZE;
		size--;
		return FIFO_SUCCESS;
	}

	uint8_t isFull() {
		return size >= FIFO_SIZE;
	}
};

struct uart2 : public sc_core::sc_module
{
    tlm_utils::simple_target_socket<uart2> tsock;
    
    interrupt_gateway *plic = 0;
    uint32_t irq_number = 0;

    FIFO rx_fifo;
    FIFO tx_fifo;

    // memory mapped configuration registers
    u_int32_t tx_data = 0;
    u_int32_t rx_data = 0;
    u_int32_t tx_ctrl = 0;
    u_int32_t rx_ctrl = 0;
    
    std::unordered_map<uint64_t, uint32_t *> addr_to_reg;

    enum {
        TX_DATA_REG_ADDR = 0x00,
        RX_DATA_REG_ADDR = 0x04,
        TX_CTRL_REG_ADDR = 0x08,
        RX_CTRL_REG_ADDR = 0x0C
    };

    SC_HAS_PROCESS(uart2);

    uart2(sc_core::sc_module_name, u_int32_t irq_number) : irq_number(irq_number)
    {
        tsock.register_b_transport(this, &uart2::transport);
        SC_THREAD(run);

        addr_to_reg = {
            {TX_DATA_REG_ADDR, &tx_data},
            {RX_DATA_REG_ADDR, &rx_data},
            {TX_CTRL_REG_ADDR, &tx_ctrl},
            {RX_CTRL_REG_ADDR, &rx_ctrl},
        };
    }

    void transport(tlm::tlm_generic_payload &trans, sc_core::sc_time &delay)
    {
        auto addr = trans.get_address();
		auto cmd = trans.get_command();
		auto len = trans.get_data_length();
		auto ptr = trans.get_data_ptr();

        assert(len == 4);
        auto it = addr_to_reg.find(addr);
        assert(it != addr_to_reg.end());    // access to non-mapped address

        if (addr <= RX_DATA_REG_ADDR && cmd == tlm::TLM_WRITE_COMMAND)
			return; // ignore write to rx_data
		
		//trigger pre read actions
        if (addr == RX_DATA_REG_ADDR && cmd == tlm::TLM_READ_COMMAND)
        {
            uint8_t c = 0;
            rx_data = (!rx_fifo.pop(&c) << 31) | c;
        }

		// actual read/write
		if (cmd == tlm::TLM_READ_COMMAND) {
			*((uint32_t *)ptr) = *it->second;
		} else if (cmd == tlm::TLM_WRITE_COMMAND) {
			*it->second = *((uint32_t *)ptr);
		} else {
			assert(false && "unsupported tlm command for uart access");
		}
    }



    void run()
    {
        while(true)
        {
            sc_core::wait(sc_core::sc_time(32000, sc_core::SC_NS));
            //create random data
            rx_fifo.push(rand() % 26 + 65);
			plic->gateway_trigger_interrupt(irq_number);
        }
    }
};


#endif //RISCV_ISA_UART2_H