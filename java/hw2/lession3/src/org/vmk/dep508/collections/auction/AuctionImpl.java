package org.vmk.dep508.collections.auction;

import java.math.BigDecimal;
import java.util.*;

//  В связи с тем, что ставки могут отменяться, логичным будет игнорировать
//  величину новой ставки(даже если она меньше предыдущей).
//  Также будем менять ставку на продукт для пользователя, который
//  уже сделал ставку(даже если она меньше предыдущей).
//  Учитывается только последняя ставка. Если пользователь сделал две ставки на один продукт,
//  то учитывается только последняя, и при удалении ставки считаем, что ставок не осталось
//  (у этого пользователя на этот продукт).
//  Продукт будет продан, если есть ставка больше установленной цены.
//  Несмотря на то, что ставка может быть снята, ценой продукта будем считать максимум
//  из всех ставок и начальной цены.
//  Так как при продаже продукта возвращается только bool, не имеет значение, кто его купил.
//  А значит, что если два пользователя делают одинаковые ставки, не имеет значения, кто
//  первый сделал ставку.

public class AuctionImpl implements Auction {
    private final Map<String, BigDecimal> products = new HashMap<>();
    private final Map<String, Map<String, BigDecimal>> bids = new HashMap<>(); // Map<product, Map<user, price>>

    @Override
    public void placeProduct(String product, BigDecimal initialPrice) {
        products.put(product, initialPrice);
    }

    @Override
    public void addBid(String user, String product, BigDecimal price) {
        if (!products.containsKey(product)) {
            throw new ProductNotFoundException("Product " + product + " not found.");
        }

        if (!bids.containsKey(product)) {
            bids.put(product, new HashMap<>());
        }

        Map<String, BigDecimal> bidForProduct = bids.get(product);
        bidForProduct.put(user, price);
    }

    @Override
    public void removeBid(String user, String product) {
        if (!bids.containsKey(product)) {
            return;
        }
        Map<String, BigDecimal> bidForProduct = bids.get(product);
        bidForProduct.remove(user);
    }

    @Override
    public boolean sellProduct(String product) {
        if (!products.containsKey(product)) {
            throw new ProductNotFoundException("Product " + product + " not found.");
        }

        BigDecimal initialPrice = products.get(product);

        if (!bids.containsKey(product)) {
            return false;
        }

        for (BigDecimal bidPrice: bids.get(product).values()) {
            if (bidPrice.compareTo(initialPrice) >= 0) {
                products.remove(product);
                bids.remove(product);
                return true;
            }
        }

        return false;
    }

    @Override
    public List<String> getProducts() {
        return new ArrayList<>(products.keySet());
    }

    @Override
    public BigDecimal getProductPrice(String product) {
        if (!products.containsKey(product)) {
            throw new ProductNotFoundException("Product " + product + " not found.");
        }
        BigDecimal price = products.get(product);
        if (!bids.containsKey(product)) {
            return price;
        }

        for (BigDecimal bidPrice: bids.get(product).values()) {
            if (bidPrice.compareTo(price) > 0) {
                price = bidPrice;
            }
        }

        return price;
    }
}
