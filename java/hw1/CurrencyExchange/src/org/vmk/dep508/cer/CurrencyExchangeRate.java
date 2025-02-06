package org.vmk.dep508.cer;

import sun.reflect.generics.reflectiveObjects.NotImplementedException;

import java.math.BigDecimal;
import java.util.Currency;

/**
 * Created by VPerov on 17.09.2018.
 */
public class CurrencyExchangeRate {
    private final Currency from;
    private final Currency to;
    private final BigDecimal rate;

    public CurrencyExchangeRate(BigDecimal rate, Currency from, Currency to) {
        this.rate = rate;
        this.from = from;
        this.to = to;
    }

    public Money convert(Money m) {
        if (!from.equals(m.getCurrency())) {
            throw new IncorrectExchangeRateException("The currency of the money object does not match the source currency.");
        }
        Currency currency = Currency.getInstance(to.getCurrencyCode());
        BigDecimal amount = rate.multiply(m.getAmount()).setScale(to.getDefaultFractionDigits(), BigDecimal.ROUND_HALF_UP);
        return new Money(currency, amount);
    }
}
