package org.vmk.dep508.cer;

import sun.reflect.generics.reflectiveObjects.NotImplementedException;

import java.math.BigDecimal;
import java.util.Currency;

public class Money {
    private Currency currency;
    private BigDecimal amount;

    public Money(Currency currency, BigDecimal amount) {
        this.currency = currency;
        this.amount = amount.setScale(this.currency.getDefaultFractionDigits());
    }

    public Currency getCurrency() {
        return currency;
    }

    public BigDecimal getAmount() {
        return amount;
    }

    public Money add(Money m) {
        if (!currency.equals(m.getCurrency())) {
            throw new DifferentCurrenciesException(
                    "Cannot add money with different currencies: "
                            + currency.getCurrencyCode() + " and "
                            + m.getCurrency().getCurrencyCode()
            );
        }
        return new Money(currency, amount.add(m.getAmount()));
    }

    public Money subtract(Money m) {
        if (!currency.equals(m.getCurrency())) {
            throw new DifferentCurrenciesException(
                    "Cannot subtract money with different currencies: "
                            + currency.getCurrencyCode() + " and "
                            + m.getCurrency().getCurrencyCode()
            );
        }
        return new Money(currency, amount.subtract(m.getAmount()));
    }

    public Money multiply(BigDecimal ratio) {
        return new Money(currency, amount.multiply(ratio).setScale(currency.getDefaultFractionDigits(), BigDecimal.ROUND_HALF_UP));
    }

    public Money devide(BigDecimal ratio) {
        return new Money(currency, amount.divide(ratio, currency.getDefaultFractionDigits(), BigDecimal.ROUND_HALF_UP));
    }
}
