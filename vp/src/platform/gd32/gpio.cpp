#include "gpio.h"

GPIO::GPIO(sc_core::sc_module_name, unsigned port) {
	tsock.register_b_transport(this, &GPIO::transport);

	router
	    .add_register_bank({
	        {GPIO_CTL0_REG_ADDR, &gpio_ctl0},
	        {GPIO_CTL1_REG_ADDR, &gpio_ctl1},
	        {GPIO_ISTAT_REG_ADDR, &gpio_istat},
	        {GPIO_OCTL_REG_ADDR, &gpio_octl},
	        {GPIO_BOP_REG_ADDR, &gpio_bop},
	        {GPIO_BC_REG_ADDR, &gpio_bc},
	        {GPIO_LOCK_REG_ADDR, &gpio_lock},
	    })
	    .register_handler(this, &GPIO::register_access_callback);

	SC_METHOD(synchronousChange);
	sensitive << asyncEvent;
	dont_initialize();

	server.setupConnection(std::to_string(port).c_str());
	server.registerOnChange(std::bind(&GPIO::asyncOnchange, this, std::placeholders::_1, std::placeholders::_2));
	serverThread = new std::thread(std::bind(&GpioServer::startAccepting, &server));
}

GPIO::~GPIO() {
	server.quit();
	if (serverThread) {
		if (serverThread->joinable()) {
			serverThread->join();
		}
		delete serverThread;
	}
}

void GPIO::asyncOnchange(gpio::PinNumber bit, gpio::Tristate val) {
	// TODO
}

void GPIO::synchronousChange() {
	// TODO
}

void GPIO::transport(tlm::tlm_generic_payload &trans, sc_core::sc_time &delay) {
	router.transport(trans, delay);
}

void GPIO::register_access_callback(const vp::map::register_access_t &r) {
	if (r.write) {
		switch (r.addr) {
			case GPIO_BOP_REG_ADDR: {
				const uint16_t set = (uint16_t)r.nv;
				const uint16_t clear = (uint16_t)(r.nv >> 16);

				uint16_t output_en = 0;
				for (int i = 0; i < 8; i++) {
					output_en |= ((gpio_ctl0 & (0b11 << (4 * i))) > 0) << i;
					output_en |= ((gpio_ctl1 & (0b11 << (4 * i))) > 0) << (i + 8);
				}

				gpio_octl &= ~(clear & output_en);
				gpio_octl |= (set & output_en);

				const uint16_t change = (set | clear) & output_en;

				for (gpio::PinNumber i = 0; i < available_pins; i++) {
					const auto bitoffs = (1l << i);
					if (bitoffs & change) {
						if (set & bitoffs)
							server.state.pins[i] = gpio::Pinstate::HIGH;
						else if (clear & bitoffs)
							server.state.pins[i] = gpio::Pinstate::LOW;
						server.pushPin(i, gpio::toTristate(server.state.pins[i]));
					}
				}
				break;
			}
			case GPIO_BC_REG_ADDR: {
				uint16_t clear = (uint16_t)r.nv;

				uint16_t output_en = 0;
				for (int i = 0; i < 8; i++) {
					output_en |= ((gpio_ctl0 & (0b11 << (4 * i))) > 0) << i;
					output_en |= ((gpio_ctl1 & (0b11 << (4 * i))) > 0) << (i + 8);
				}

				gpio_octl &= ~(clear & output_en);

				for (gpio::PinNumber i = 0; i < available_pins; i++) {
					const auto bitoffs = (1l << i);
					if (bitoffs & clear & output_en) {
						server.state.pins[i] = gpio::Pinstate::LOW;
						server.pushPin(i, gpio::toTristate(server.state.pins[i]));
					}
				}
				break;
			}
			case GPIO_OCTL_REG_ADDR: {
				for (gpio::PinNumber i = 0; i < available_pins; i++) {
					const auto bitoffs = (1l << i);
					server.state.pins[i] = (uint16_t)r.nv & bitoffs ? gpio::Pinstate::HIGH : gpio::Pinstate::LOW;
					server.pushPin(i, gpio::toTristate(server.state.pins[i]));
				}
				break;
			}

			default:
				break;
		}
	}
	r.fn();
}
